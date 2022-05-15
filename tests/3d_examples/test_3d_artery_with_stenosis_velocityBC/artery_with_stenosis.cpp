/* --------------------------------------------------------------------------- *
 *                       SPHinXsys: 3D artery with stenosis example                        *
 * ----------------------------------------------------------------------------*
 * This is the test case for  3D artery with stenosis example by velocity boundary condition.  *
 * ---------------------------------------------------------------------------*/
/**
 * @brief 	SPHinXsys Library.
 */
#include "sphinxsys.h"

using namespace SPH;

// for geometry
Real resolution_ref = 0.0003; // particle spacing

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(-0.1, -0.002, -0.014), Vecd(0.002, 0.003, -.007));
// for material properties of the fluid
Real rho0_f = 1050.0; // kg/m3
Real U_f = 0.02;	  // m/s
Real U_max = 0.1;
Real L_f = 3e-3; // reference length
Real c_f = 10.0 * U_max;
Real mu_f = 3.5e-3; // in Pa-s

// boundary faces
Vecd inlet_center(-0.077000, 0.00055226, -0.010811);
Vecd inlet_direction(1, 0, 0);
StdVec<Vecd> inlet_points{
	inlet_center + Vecd(0, 0.002, 0),
	inlet_center - Vecd(0, 0.002, 0),
	inlet_center + Vecd(0, 0, 0.002),
	inlet_center - Vecd(0, 0, 0.002)};
SegmentFace inlet_face(inlet_points, inlet_direction, inlet_center);

Vecd outlet_center(-0.022000, 0.00055226, -0.010811);
Vecd outlet_direction(-1, 0, 0);
StdVec<Vecd> outlet_points{
	outlet_center + Vecd(0, 0.002, 0),
	outlet_center - Vecd(0, 0.002, 0),
	outlet_center + Vecd(0, 0, 0.002),
	outlet_center - Vecd(0, 0, 0.002)};
SegmentFace outlet_face(outlet_points, outlet_direction, outlet_center);

//	import the fluid body
class Lumen : public ComplexShape
{
public:
	explicit Lumen(const std::string &shape_name) : ComplexShape(shape_name)
	{
		std::string file_name = "./input/fluid.stl";		  // in millimeter
		add<TriangleMeshShapeSTL>(file_name, Vecd(0), 0.001); // millimeter to meter
	}
};

//	import the static solid wall boundary
class Artery : public ComplexShape
{
public:
	explicit Artery(const std::string &shape_name) : ComplexShape(shape_name)
	{
		std::string file_name = "./input/solid.stl";
		add<TriangleMeshShapeSTL>(file_name, Vecd(0), 0.001);
	}
};

/**
 * @brief define velocity buffer
 */
class EmitterBufferInflowConditionWithFace : public VelocityInflowConditionWithFace
{
	Real u_ave_, u_ref_, t_ref_;

public:
	EmitterBufferInflowConditionWithFace(FluidBody &body, BodyRegionByCellsWithFace &body_part)
		: VelocityInflowConditionWithFace(body, body_part),
		  u_ave_(0), u_ref_(U_f), t_ref_(2.0) {}

	Vecd defineVelocityProfile(Vecd &position, Vecd &velocity) override
	{
		Vecd vel = u_ave_ * inlet_direction.normalize();
		return vel;
	}

	void setupDynamics(Real dt = 0.0) override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		u_ave_ = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
	}
};

// the main program
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;
	GlobalStaticVariables::physical_time_ = 0;

	/* relax wall particles at first */
	system.run_particle_relaxation_ = true;
	system.reload_particles_ = false;

	// handle command line arguments
	system.handleCommandlineOptions(ac, av);
	InOutput in_output(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FluidBody lumen(system, makeShared<Lumen>("Lumen"));
	lumen.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	lumen.generateParticles<ParticleGeneratorLattice>();

	SolidBody artery(system, makeShared<Artery>("Artery"));
	artery.defineBodyLevelSetShape()->writeLevelSet(artery);
	artery.defineParticlesAndMaterial<SolidParticles, Solid>();
	(!system.run_particle_relaxation_ && system.reload_particles_)
		? artery.generateParticles<ParticleGeneratorReload>(in_output, artery.getBodyName())
		: artery.generateParticles<ParticleGeneratorLattice>();
	SimpleDynamics<NormalDirectionFromBodyShape> artery_normal_direction(artery);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner lumen_inner(lumen);
	ComplexBodyRelation lumen_complex_relation(lumen_inner, {&artery});
	//----------------------------------------------------------------------
	// Relaxing solid particles.
	//----------------------------------------------------------------------
	if (system.run_particle_relaxation_)
	{
		/* body topology only for particle relation */
		BodyRelationInner wall_bound_model_inner(artery);
		/* Random reset the insert body particle position */
		RandomizePartilePosition random_wall_particles(artery);
		/* write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_inserted_body_to_vtp(in_output, {&artery});
		/* write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(in_output, {&artery});
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(wall_bound_model_inner, true);
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_wall_particles.parallel_exec(0.25);
		relaxation_step_inner.surface_bounding_.parallel_exec();
		write_inserted_body_to_vtp.writeToFile(0);
		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner.parallel_exec();
			ite_p++;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the imported model N = " << ite_p << "\n";
				write_inserted_body_to_vtp.writeToFile(ite_p);
			}
		}
		cout << "The physics relaxation process of artery body finish ! \n";

		/* output results. */
		write_particle_reload_files.writeToFile(0);
		return 0;
	}
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	/** Initialize particle acceleration. */
	TimeStepInitialization initialize_a_fluid_step(lumen);
	/* Define particle emitter at inlet */
	BodyRegionByParticleWithFace emitter(lumen, inlet_face, 4);
	InflowInjectingWithFace emitter_inflow_injecting(lumen, emitter, 20);

	/** Define inflow condition. */
	BodyRegionByCellsWithFace emitter_buffer(lumen, inlet_face, 4);
	EmitterBufferInflowConditionWithFace emitter_buffer_inflow_condition(lumen, emitter_buffer);

	/** time-space method to detect surface particles. */
	fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex
		inlet_outlet_surface_particle_indicator(lumen_complex_relation);
	/** Evaluation of density by freestream approach. */
	fluid_dynamics::DensitySummationFreeStreamComplex update_density_by_summation(lumen_complex_relation);
	/** We can output a method-specific particle data for debug */
	lumen.addBodyStateForRecording<Real>("Pressure");
	lumen.addBodyStateForRecording<int>("SurfaceIndicator");
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize get_fluid_advection_time_step_size(lumen, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(lumen);
	/** Pressure relaxation. */
	fluid_dynamics::PressureRelaxationWithWall pressure_relaxation(lumen_complex_relation);
	/** Density relaxation. */
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation(lumen_complex_relation);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(lumen_complex_relation);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityCorrectionComplex transport_velocity_correction(lumen_complex_relation);
	/** Define outlet face to */
	BodyRegionByCellsWithFace face1(lumen, outlet_face);
	/** recycle real fluid particle to buffer particles at outlet. */
	DoNothingConditionWithFace outflow1(lumen, face1);

	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_body_states(in_output, system.real_bodies_);
	RestartIO restart_io(in_output, system.real_bodies_);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	artery_normal_direction.parallel_exec();
	inlet_outlet_surface_particle_indicator.parallel_exec();
	//----------------------------------------------------------------------
	//	Load restart file if necessary.
	//----------------------------------------------------------------------
	/** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.restart_step_);
		lumen.updateCellLinkedList();
		lumen_complex_relation.updateConfiguration();
	}
	//----------------------------------------------------------------------
	//	Setup computing and initial conditions.
	//----------------------------------------------------------------------
	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 1.0;		   /**< End time. */
	Real D_Time = End_Time / 50.0; /**< Time stamps for output of body states. */
	Real dt = 0.0;				   /**< Default acoustic time step sizes. */
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_body_states.writeToFile();
	//----------------------------------------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			initialize_a_fluid_step.parallel_exec();
			Real Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_density_by_summation.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec(Dt);

			/** Dynamics including pressure relaxation. */
			Real relaxation_time = 0.0;
			int iii = 0;
			while (relaxation_time < Dt)
			{
				dt = SMIN(get_fluid_time_step_size.parallel_exec(), Dt - relaxation_time);
				pressure_relaxation.parallel_exec(dt);
				emitter_buffer_inflow_condition.parallel_exec();
				density_relaxation.parallel_exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				++iii;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	dt = " << dt << "  Dt/dt= " << iii << "\n";

				if (number_of_iterations % restart_output_interval == 0 && number_of_iterations != system.restart_step_)
					restart_io.writeToFile(Real(number_of_iterations));
			}
			number_of_iterations++;

			/** inflow injecting */
			emitter_inflow_injecting.exec();

			/* check particle at outlet */
			outflow1.parallel_exec();

			/** Update cell linked list and configuration. */
			lumen.updateCellLinkedList();
			lumen_complex_relation.updateConfiguration();
			inlet_outlet_surface_particle_indicator.parallel_exec();
		}

		tick_count t2 = tick_count::now();
		write_body_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
			  << " seconds." << std::endl;

	return 0;
}
