/**
 * @file 	windows_frame_diffusion_D2.h
 * @brief 	Head file of windows_frame_diffusion_D2.cpp.
 * @author	Haotian Ji, Dong Wu, Chi Zhang and Xiangyu Hu
 */
#ifndef WINDOWS_FRAME_DIFFUSION_D2_H
#define WINDOWS_FRAME_DIFFUSION_D2_H

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 0.3; // Unit in m.
Real H = 0.129;
Real resolution_ref = 0.001;
Real BW = resolution_ref * 2.0;
Real particle_volume = pow(resolution_ref, 2) * 1; // consider unit length in z direction for 2d case
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(L + BW, H + BW));
//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
// Const for calculating air cavities conductivity
Real C1 = 0.025;        // Unit W/(m*K)
Real C3 = 1.57;         // Unit W/(m2*K)
Real C4 = 2.11;         // Unit W/(m*K)
Real alum_cond = 160;   // Unit W/(m*k), aluminium conductivity
Real wood_cond = 0.13;   // wood conductivity
Real epdm_cond = 0.25;  // epdm conductivity
Real pane_cond = 0.035; // insulation panel conductivity

Real getACConductivity(Real b, Real d, Real A) // calculate air cavities conductivity
{
    Real b_equal = sqrt(A * b / d);
    Real d_equal = sqrt(A * d / b);

    Real ha = 0.0;
    if (b_equal < 0.005)
    {
        ha = C1 / d_equal;
    }
    else
    {
        ha = SMAX(C1 / d_equal, C3);
    }

    Real hr = C4 * (1 - d_equal / b_equal + sqrt(1 + pow(d_equal / b_equal, 2)));
    Real Rs = 1 / (ha + hr);
    Real cond = d_equal / Rs;

    return cond;
}

/*---Geometric parameter of air cavities, unit m and m2, d is the cavity dimension
in the heat flow rate direction, b is the cavity dimension perpendicular to the heat
flow rate direction. The Area parameters for non-rectangular are given, for rectangular
formulas are given.---*/
Real d1 = 0.006;Real b1 = 0.076;Real A1 = d1 * b1;
Real d2 = 0.022;Real b2 = 0.010;Real A2 = d2 * b2;
Real d3 = 0.027;Real b3 = 0.030;Real A3 = 0.000589;
Real d4 = 0.015;Real b4 = 0.013;Real A4 = 0.000181;
Real d5 = 0.026;Real b5 = 0.005;Real A5 = d5 * b5;
Real d6 = 0.077;Real b6 = 0.051;Real A6 = 0.001047;

// unventilated air cavities conductivity
Real ac1_cond = getACConductivity(b1, d1, A1);
Real ac2_cond = getACConductivity(b2, d2, A2);
Real ac3_cond = getACConductivity(b3, d3, A3);
Real ac4_cond = getACConductivity(b4, d4, A4);
Real ac5_cond = getACConductivity(b5, d5, A5);
Real ac6_cond = getACConductivity(b6, d6, A6);
//----------------------------------------------------------------------
//	Initial and boundary conditions.
//----------------------------------------------------------------------
// Temperature initialization,unit Celsius
Real initial_temperature = 10.0;
Real T_infinity_e = 0.0;
Real T_infinity_i = 20.0;

//Surface resistance, unit (K*m2)/W
Real rs_i = 0.13;
Real rs_i_increased = 0.20;
Real rs_e = 0.04;

// Convection coefficient, unit W/(K*m2)
Real convection_e = 1 / rs_e;
Real convection_i = 1 / rs_i;
Real convection_i_decreased = 1 / rs_i_increased;
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
MultiPolygon createOverallStructureBody()
{
    std::vector<Vecd> overallStructureDomainShape;
    overallStructureDomainShape.push_back(Vecd(0.0, 0.005));
    overallStructureDomainShape.push_back(Vecd(0.0, 0.105));
    overallStructureDomainShape.push_back(Vecd(0.047, 0.105));
    overallStructureDomainShape.push_back(Vecd(0.047, 0.124));
    overallStructureDomainShape.push_back(Vecd(0.11, 0.124));
    overallStructureDomainShape.push_back(Vecd(0.11, 0.105));
    overallStructureDomainShape.push_back(Vecd(0.3, 0.105));
    overallStructureDomainShape.push_back(Vecd(0.3, 0.047));
    overallStructureDomainShape.push_back(Vecd(0.11, 0.047));
    overallStructureDomainShape.push_back(Vecd(0.11, 0.033));
    overallStructureDomainShape.push_back(Vecd(0.09, 0.033));
    overallStructureDomainShape.push_back(Vecd(0.09, 0.005));
    overallStructureDomainShape.push_back(Vecd(0.0, 0.005));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(overallStructureDomainShape, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createInternalAirBody()
{
    std::vector<Vecd> internalAirDomainShape;
    internalAirDomainShape.push_back(Vecd(0.0, 0.105));
    internalAirDomainShape.push_back(Vecd(0.0, 0.11));
    internalAirDomainShape.push_back(Vecd(0.028, 0.11));
    internalAirDomainShape.push_back(Vecd(0.028, 0.129));
    internalAirDomainShape.push_back(Vecd(0.129, 0.129));
    internalAirDomainShape.push_back(Vecd(0.129, 0.11));
    internalAirDomainShape.push_back(Vecd(0.3, 0.11));
    internalAirDomainShape.push_back(Vecd(0.3, 0.105));
    internalAirDomainShape.push_back(Vecd(0.11, 0.105));
    internalAirDomainShape.push_back(Vecd(0.11, 0.124));
    internalAirDomainShape.push_back(Vecd(0.047, 0.124));
    internalAirDomainShape.push_back(Vecd(0.047, 0.105));
    internalAirDomainShape.push_back(Vecd(0.0, 0.105));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(internalAirDomainShape, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createDecreasedInternalConvectionBody()
{
    std::vector<Vecd> decreasedInternalConvectionDomainShape1;
    decreasedInternalConvectionDomainShape1.push_back(Vecd(0.028, 0.105));
    decreasedInternalConvectionDomainShape1.push_back(Vecd(0.047, 0.124));
    decreasedInternalConvectionDomainShape1.push_back(Vecd(0.047, 0.105));
    decreasedInternalConvectionDomainShape1.push_back(Vecd(0.028, 0.105));

    std::vector<Vecd> decreasedInternalConvectionDomainShape2;
    decreasedInternalConvectionDomainShape2.push_back(Vecd(0.110, 0.105));
    decreasedInternalConvectionDomainShape2.push_back(Vecd(0.110, 0.124));
    decreasedInternalConvectionDomainShape2.push_back(Vecd(0.129, 0.105));
    decreasedInternalConvectionDomainShape2.push_back(Vecd(0.110, 0.105));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(decreasedInternalConvectionDomainShape1, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(decreasedInternalConvectionDomainShape2, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createExternalAirBody()
{
    std::vector<Vecd> externalAirDomainShape;
    externalAirDomainShape.push_back(Vecd(0.0, 0.0));
    externalAirDomainShape.push_back(Vecd(0.0, 0.005));
    externalAirDomainShape.push_back(Vecd(0.09, 0.005));
    externalAirDomainShape.push_back(Vecd(0.09, 0.033));
    externalAirDomainShape.push_back(Vecd(0.11, 0.033));
    externalAirDomainShape.push_back(Vecd(0.11, 0.047));
    externalAirDomainShape.push_back(Vecd(0.3, 0.047));
    externalAirDomainShape.push_back(Vecd(0.3, 0.042));
    externalAirDomainShape.push_back(Vecd(0.115, 0.042));
    externalAirDomainShape.push_back(Vecd(0.115, 0.028));
    externalAirDomainShape.push_back(Vecd(0.095, 0.028));
    externalAirDomainShape.push_back(Vecd(0.095, 0.0));
    externalAirDomainShape.push_back(Vecd(0.0, 0.0));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(externalAirDomainShape, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createWoodBody()
{
    std::vector<Vecd> woodDomainShape1;
    woodDomainShape1.push_back(Vecd(0.0, 0.013));
    woodDomainShape1.push_back(Vecd(0.0, 0.105));
    woodDomainShape1.push_back(Vecd(0.059, 0.105));
    woodDomainShape1.push_back(Vecd(0.059, 0.031));
    woodDomainShape1.push_back(Vecd(0.076, 0.031));
    woodDomainShape1.push_back(Vecd(0.076, 0.013));
    woodDomainShape1.push_back(Vecd(0.0, 0.013));

    std::vector<Vecd> woodDomainShape2;
    woodDomainShape2.push_back(Vecd(0.047, 0.108));
    woodDomainShape2.push_back(Vecd(0.047, 0.124));
    woodDomainShape2.push_back(Vecd(0.11, 0.124));
    woodDomainShape2.push_back(Vecd(0.11, 0.108));
    woodDomainShape2.push_back(Vecd(0.09, 0.108));
    woodDomainShape2.push_back(Vecd(0.09, 0.082));
    woodDomainShape2.push_back(Vecd(0.11, 0.082));
    woodDomainShape2.push_back(Vecd(0.11, 0.07));
    woodDomainShape2.push_back(Vecd(0.063, 0.07));
    woodDomainShape2.push_back(Vecd(0.063, 0.108));
    woodDomainShape2.push_back(Vecd(0.047, 0.108));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(woodDomainShape1, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(woodDomainShape2, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createEPDMBody()
{
    std::vector<Vecd> epdmDomainShape1;
    epdmDomainShape1.push_back(Vecd(0.047, 0.105));
    epdmDomainShape1.push_back(Vecd(0.047, 0.108));
    epdmDomainShape1.push_back(Vecd(0.059, 0.108));
    epdmDomainShape1.push_back(Vecd(0.059, 0.105));
    epdmDomainShape1.push_back(Vecd(0.047, 0.105));

    std::vector<Vecd> epdmDomainShape2;
    epdmDomainShape2.push_back(Vecd(0.095, 0.105));
    epdmDomainShape2.push_back(Vecd(0.095, 0.108));
    epdmDomainShape2.push_back(Vecd(0.11, 0.108));
    epdmDomainShape2.push_back(Vecd(0.11, 0.105));
    epdmDomainShape2.push_back(Vecd(0.095, 0.105));

    std::vector<Vecd> epdmDomainShape3;
    epdmDomainShape3.push_back(Vecd(0.095, 0.082));
    epdmDomainShape3.push_back(Vecd(0.095, 0.085));
    epdmDomainShape3.push_back(Vecd(0.11, 0.085));
    epdmDomainShape3.push_back(Vecd(0.11, 0.082));
    epdmDomainShape3.push_back(Vecd(0.095, 0.082));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(epdmDomainShape1, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(epdmDomainShape2, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(epdmDomainShape3, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createPanelBody()
{
    std::vector<Vecd> panelDomainShape;
    panelDomainShape.push_back(Vecd(0.095, 0.085));
    panelDomainShape.push_back(Vecd(0.095, 0.105));
    panelDomainShape.push_back(Vecd(0.3, 0.105));
    panelDomainShape.push_back(Vecd(0.3, 0.047));
    panelDomainShape.push_back(Vecd(0.11, 0.047));
    panelDomainShape.push_back(Vecd(0.11, 0.085));
    panelDomainShape.push_back(Vecd(0.095, 0.085));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(panelDomainShape, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createACBody1()
{
    std::vector<Vecd> acDomainShape1;
    acDomainShape1.push_back(Vecd(0.0, 0.007));
    acDomainShape1.push_back(Vecd(0.000, 0.013));
    acDomainShape1.push_back(Vecd(0.076, 0.013));
    acDomainShape1.push_back(Vecd(0.076, 0.007));
    acDomainShape1.push_back(Vecd(0.0, 0.007));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(acDomainShape1, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createACBody2()
{
    std::vector<Vecd> acDomainShape2;
    acDomainShape2.push_back(Vecd(0.078, 0.007));
    acDomainShape2.push_back(Vecd(0.078, 0.029));
    acDomainShape2.push_back(Vecd(0.088, 0.029));
    acDomainShape2.push_back(Vecd(0.088, 0.007));
    acDomainShape2.push_back(Vecd(0.078, 0.007));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(acDomainShape2, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createACBody3()
{
    std::vector<Vecd> acDomainShape3;
    acDomainShape3.push_back(Vecd(0.078, 0.035));
    acDomainShape3.push_back(Vecd(0.078, 0.062));
    acDomainShape3.push_back(Vecd(0.095, 0.062));
    acDomainShape3.push_back(Vecd(0.095, 0.045));
    acDomainShape3.push_back(Vecd(0.108, 0.045));
    acDomainShape3.push_back(Vecd(0.108, 0.035));
    acDomainShape3.push_back(Vecd(0.078, 0.035));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(acDomainShape3, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createACBody4()
{
    std::vector<Vecd> acDomainShape4;
    acDomainShape4.push_back(Vecd(0.097, 0.047));
    acDomainShape4.push_back(Vecd(0.097, 0.062));
    acDomainShape4.push_back(Vecd(0.108, 0.062));
    acDomainShape4.push_back(Vecd(0.108, 0.055));
    acDomainShape4.push_back(Vecd(0.11, 0.055));
    acDomainShape4.push_back(Vecd(0.11, 0.047));
    acDomainShape4.push_back(Vecd(0.097, 0.047));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(acDomainShape4, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createACBody5()
{
    std::vector<Vecd> acDomainShape5;
    acDomainShape5.push_back(Vecd(0.09, 0.082));
    acDomainShape5.push_back(Vecd(0.09, 0.108));
    acDomainShape5.push_back(Vecd(0.095, 0.108));
    acDomainShape5.push_back(Vecd(0.095, 0.082));
    acDomainShape5.push_back(Vecd(0.09, 0.082));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(acDomainShape5, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createACBody6()
{
    std::vector<Vecd> acDomainShape6;
    acDomainShape6.push_back(Vecd(0.059, 0.031));
    acDomainShape6.push_back(Vecd(0.059, 0.108));
    acDomainShape6.push_back(Vecd(0.063, 0.108));
    acDomainShape6.push_back(Vecd(0.063, 0.07));
    acDomainShape6.push_back(Vecd(0.11, 0.07));
    acDomainShape6.push_back(Vecd(0.11, 0.064));
    acDomainShape6.push_back(Vecd(0.076, 0.064));
    acDomainShape6.push_back(Vecd(0.076, 0.033));
    acDomainShape6.push_back(Vecd(0.09, 0.033));
    acDomainShape6.push_back(Vecd(0.09, 0.031));
    acDomainShape6.push_back(Vecd(0.059, 0.031));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(acDomainShape6, ShapeBooleanOps::add);

    return multi_polygon;
}
//----------------------------------------------------------------------
//	Setup diffusion material properties.
//----------------------------------------------------------------------
class DiffusionMaterial : public DiffusionReaction<Solid>
{
  public:
    DiffusionMaterial() : DiffusionReaction<Solid>({"Phi"}, SharedPtr<NoReaction>())
    {
        initializeAnDiffusion<LocalIsotropicDiffusion>("Phi", "Phi", alum_cond);
    }
};

class WallMaterial : public DiffusionReaction<Solid>
{
  public:
    WallMaterial() : DiffusionReaction<Solid>({"Phi"}, SharedPtr<NoReaction>())
    {
        initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi", 110);
        // We should use aluminium conductivity to define the time step in principle cause it's the largest,
        // 110 comes from iteration, which can accelarate the simulation and also remain stability.
    }
};
using DiffusionParticles = DiffusionReactionParticles<SolidParticles, DiffusionMaterial>;
using WallParticles = DiffusionReactionParticles<SolidParticles, WallMaterial>;
//----------------------------------------------------------------------
//	Application dependent initial condition.
//----------------------------------------------------------------------
template <class DynamicsIdentifier, class ParticlesType>
class LocalQuantityDefinition
    : public BaseLocalDynamics<DynamicsIdentifier>,
      public DiffusionReactionSimpleData<ParticlesType>
{
  public:
    LocalQuantityDefinition(DynamicsIdentifier &identifier)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier),
          DiffusionReactionSimpleData<ParticlesType>(identifier.getSPHBody()){};
    virtual ~LocalQuantityDefinition(){};
};
class ThermalConductivityInitialization
    : public LocalQuantityDefinition<BodyPartByParticle, DiffusionParticles>
{
  protected:
    StdLargeVec<Real> &thermal_conductivity;
    Real local_diff;

  public:
    explicit ThermalConductivityInitialization(BodyPartByParticle &body_part, Real local_diff)
        : LocalQuantityDefinition<BodyPartByParticle, DiffusionParticles>(body_part),
          thermal_conductivity(*(particles_->getVariableByName<Real>("ThermalConductivity"))),
          local_diff(local_diff){};

    void update(size_t index_i, Real dt)
    {
        thermal_conductivity[index_i] = local_diff;
    }
};
class LocalConvectionInitialization
    : public LocalQuantityDefinition<BodyPartByParticle, WallParticles>
{
  protected:
    StdLargeVec<Real> &convection_;
    Real local_convection;

  public:
    explicit LocalConvectionInitialization(BodyPartByParticle &body_part, Real local_convection)
        : LocalQuantityDefinition<BodyPartByParticle, WallParticles>(body_part),
          convection_(*(this->particles_->template getVariableByName<Real>("Convection"))),
          local_convection(local_convection){};

    void update(size_t index_i, Real dt)
    {
        convection_[index_i] = local_convection;
    }
};
class LocalHeatTransferConvection
    : public LocalQuantityDefinition<BodyPartByParticle, WallParticles>
{
  protected:
    StdLargeVec<Real> &ht_convection_;
    Real local_ht_convection;

  public:
    explicit LocalHeatTransferConvection(BodyPartByParticle &body_part, Real local_ht_convection)
        : LocalQuantityDefinition<BodyPartByParticle, WallParticles>(body_part),
          ht_convection_(*(this->particles_->template getVariableByName<Real>("HT_Convection"))),
          local_ht_convection(local_ht_convection){};

    void update(size_t index_i, Real dt)
    {
        ht_convection_[index_i] = local_ht_convection;
    }
};
class DiffusionInitialCondition
    : public DiffusionReactionInitialCondition<DiffusionParticles>
{
  protected:
    size_t phi_;

  public:
    explicit DiffusionInitialCondition(SPHBody &sph_body)
        : DiffusionReactionInitialCondition<DiffusionParticles>(sph_body)
    {
        phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
    };

    void update(size_t index_i, Real dt)
    {
        all_species_[phi_][index_i] = initial_temperature;
    };
};
class RobinWallBoundaryInitialCondition
    : public DiffusionReactionInitialCondition<WallParticles>
{
  protected:
    size_t phi_;
    StdLargeVec<Real> &convection_;
    Real &T_infinity_;

  public:
    explicit RobinWallBoundaryInitialCondition(SolidBody &diffusion_body)
        : DiffusionReactionInitialCondition<WallParticles>(diffusion_body),
          convection_(*(this->particles_->template getVariableByName<Real>("Convection"))),
          T_infinity_(*(this->particles_->template getGlobalVariableByName<Real>("T_infinity")))
    {
        phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
    }

    void update(size_t index_i, Real dt)
    {
        all_species_[phi_][index_i] = -0.0;

        if (pos_[index_i][1] >= 0.105)
        {
            convection_[index_i] = convection_i;
            T_infinity_ = T_infinity_i;
        }
        if (pos_[index_i][1] <= 0.047)
        {
            convection_[index_i] = convection_e;
            T_infinity_ = T_infinity_e;
        }
    }
};
class ExternalHeatTransferInitialCondition
    : public DiffusionReactionInitialCondition<WallParticles>
{
  protected:
    size_t phi_;
    StdLargeVec<Real> &ht_convection_;
    Real &ht_T_infinity_;

  public:
    explicit ExternalHeatTransferInitialCondition(SolidBody &diffusion_body)
        : DiffusionReactionInitialCondition<WallParticles>(diffusion_body),
          ht_convection_(*(this->particles_->template getVariableByName<Real>("HT_Convection"))),
          ht_T_infinity_(*(this->particles_->template getGlobalVariableByName<Real>("HT_T_infinity")))
    {
        phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
    }

    void update(size_t index_i, Real dt)
    {
        ht_convection_[index_i] = convection_e;
        ht_T_infinity_ = T_infinity_e;
    }
};
class InternalHeatTransferInitialCondition
    : public DiffusionReactionInitialCondition<WallParticles>
{
  protected:
    size_t phi_;
    StdLargeVec<Real> &ht_convection_;
    Real &ht_T_infinity_;

  public:
    explicit InternalHeatTransferInitialCondition(SolidBody &diffusion_body)
        : DiffusionReactionInitialCondition<WallParticles>(diffusion_body),
          ht_convection_(*(this->particles_->template getVariableByName<Real>("HT_Convection"))),
          ht_T_infinity_(*(this->particles_->template getGlobalVariableByName<Real>("HT_T_infinity")))
    {
        phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
    }

    void update(size_t index_i, Real dt)
    {
        ht_convection_[index_i] = convection_i;
        ht_T_infinity_ = T_infinity_i;
    }
};

template <typename... ControlTypes>
class RobinFlux; /*Calculate heat transfer flux of Robin boundary condition*/
template <class ParticlesType, class ContactParticlesType, class ContactKernelGradientType>
class DiffusionRelaxation<RobinFlux<ParticlesType, ContactParticlesType, ContactKernelGradientType>>
    : public DiffusionRelaxation<BaseContact, ParticlesType, ContactParticlesType, ContactKernelGradientType>
{
    StdLargeVec<Vecd> &n_;
    StdVec<StdLargeVec<Vecd> *> ht_n_;
    StdVec<Real *> ht_T_infinity_; 
    StdVec<StdLargeVec<Real> *> ht_flux_;
    StdVec<StdLargeVec<Real> *> ht_convection_; // ht for heat transfer

  public:
    explicit DiffusionRelaxation(BaseContactRelation &contact_relation)
        : DiffusionRelaxation<BaseContact, ParticlesType, ContactParticlesType, ContactKernelGradientType>(contact_relation),
          n_(this->particles_->n_)
    {
       ht_flux_.resize(this->all_diffusions_.size());
       ht_T_infinity_.resize(this->all_diffusions_.size());
       ht_convection_.resize(this->contact_particles_.size());

       for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
       {
            ht_flux_[m] = this->particles_->template registerSharedVariable<Real>("HT_Flux");
            
            for (size_t k = 0; k != this->contact_particles_.size(); ++k)
            {
                ht_n_.push_back(&(this->contact_particles_[k]->n_));
                ht_convection_[k] = this->contact_particles_[k]->template registerSharedVariable<Real>("HT_Convection");
                ht_T_infinity_[m] = this->contact_particles_[k]->template registerGlobalVariable<Real>("HT_T_infinity");
            }
       }
    };
    virtual ~DiffusionRelaxation(){};

    void update(size_t index_i, Real dt)
    {
        for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
        {
            for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
            {
                Real ht_flux_ = 0.0;

                StdLargeVec<Vecd> &n_k = *(ht_n_[k]);
                StdLargeVec<Real> &ht_convection_k = *(ht_convection_[k]);
                Real &ht_T_infinity_k = *(ht_T_infinity_[k]);

                Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
                for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                {
                    size_t index_j = contact_neighborhood.j_[n];
                    Real dW_ijV_j_ = contact_neighborhood.dW_ijV_j_[n];
                    Vecd &e_ij = contact_neighborhood.e_ij_[n];

                    const Vecd &grad_ijV_j = this->contact_kernel_gradients_[k](index_i, index_j, dW_ijV_j_, e_ij);
                    Vecd n_ij = n_[index_i] - n_k[index_j];
                    Real area_ij_Robin = grad_ijV_j.dot(n_ij);

                    Real phi_ij = ht_T_infinity_k - (*this->diffusion_species_[m])[index_i];
                    ht_flux_ += ht_convection_k[index_j] * phi_ij * area_ij_Robin * particle_volume;
                }
                (*this->ht_flux_[m])[index_i] = ht_flux_;
            }
        }
    };
};

using DiffusionBodyRelaxation = DiffusionBodyRelaxationComplex<
    DiffusionParticles, WallParticles, KernelGradientInner, KernelGradientContact, Robin, Robin>;
using RobinFluxCalculation = SimpleDynamics<DiffusionRelaxation<
    RobinFlux<DiffusionParticles, WallParticles, KernelGradientContact>>>;

#endif // WINDOWS_FRAME_DIFFUSION_D2_H