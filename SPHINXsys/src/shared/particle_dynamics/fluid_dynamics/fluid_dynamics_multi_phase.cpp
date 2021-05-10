/**
 * @file 	fluid_dynamics_multi_phase.cpp
 * @author	Chi ZHang and Xiangyu Hu
 */

#include "fluid_dynamics_multi_phase.h"

//=================================================================================================//
namespace SPH
{
//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		ViscousAccelerationMultiPhase::ViscousAccelerationMultiPhase(BaseInnerBodyRelation* inner_relation,
			BaseContactBodyRelation* contact_relation)
		: ViscousAccelerationInner(inner_relation), MultiPhaseContactData(contact_relation)
		{
			if (inner_relation->sph_body_ != contact_relation->sph_body_)
			{
				std::cout << "\n Error: the two body_realtions do not have the same source body!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}

			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
				contact_vel_n_.push_back(&(contact_particles_[k]->vel_n_));
			}
		}
		//=================================================================================================//
        ViscousAccelerationMultiPhase::ViscousAccelerationMultiPhase(ComplexBodyRelation* complex_relation) 
		: ViscousAccelerationMultiPhase(complex_relation->inner_relation_, complex_relation->contact_relation_) {}
		//=================================================================================================//
		void ViscousAccelerationMultiPhase::Interaction(size_t index_i, Real dt)
		{
			ViscousAccelerationInner::Interaction(index_i, dt);

			Real rho_i = this->rho_n_[index_i];
			Vecd& vel_i = this->vel_n_[index_i];

			Vecd acceleration(0), vel_derivative(0);
			for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
			{
				Real mu_j = this->contact_material_[k]->ReferenceViscosity();
				StdLargeVec<Real>& Vol_k = *(this->contact_Vol_[k]);
				StdLargeVec<Vecd>& vel_k = *(this->contact_vel_n_[k]);
				Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];

					vel_derivative = 2.0*(vel_i - vel_k[index_j]) / (r_ij + 0.01 * this->smoothing_length_);
					Real mu_ij = 2.0 * this->mu_ * mu_j / (this->mu_ + mu_j);
					acceleration += 2.0 * mu_ij * vel_derivative 
								  * contact_neighborhood.dW_ij_[n] * Vol_k[index_j] / rho_i;
				}
			}

			dvel_dt_others_[index_i] += acceleration;
		}
		//=================================================================================================//
		MultiPhaseColorFunctionGradient::
			MultiPhaseColorFunctionGradient(BaseContactBodyRelation* contact_relation) : 
			InteractionDynamics(contact_relation->sph_body_), MultiPhaseData(contact_relation),
			rho_0_(particles_->rho_0_), Vol_(particles_->Vol_),
			pos_div_(*particles_->getVariableByName<indexScalar, Real>("PositionDivergence")), 
			surface_indicator_(particles_->surface_indicator_),
			color_grad_(*particles_->createAVariable<indexVector, Vecd>("ColorGradient")),
			surface_norm_(*particles_->createAVariable<indexVector, Vecd>("SurfaceNormal"))
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k) 
			{	Real rho_0_k = contact_particles_[k]->rho_0_;
				contact_rho_0_.push_back(rho_0_k);
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
			}
		}
		//=================================================================================================//
		void MultiPhaseColorFunctionGradient::Interaction(size_t index_i, Real dt)
		{
			Real vol_i = Vol_[index_i];
			Real pos_div = 0.0;
			Vecd gradient(0.0);
			if(surface_indicator_[index_i])
			{
				for (size_t k = 0; k < contact_configuration_.size(); ++k)
				{
					Real rho_0_k = contact_rho_0_[k];
					StdLargeVec<Real>& contact_vol_k = *(contact_Vol_[k]);
					Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
					for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
					{
						size_t index_j = contact_neighborhood.j_[n];
						pos_div -= contact_neighborhood.dW_ij_[n] * contact_neighborhood.r_ij_[n] * contact_vol_k[index_j];
						/** Norm of interface.*/
						Real rho_ij = rho_0_ / (rho_0_ + rho_0_k);
						Real area_ij = (vol_i * vol_i + contact_vol_k[index_j] * contact_vol_k[index_j]) * contact_neighborhood.dW_ij_[n];
						gradient += rho_ij * area_ij * contact_neighborhood.e_ij_[n] / vol_i;
					}
				}
			}
			color_grad_[index_i] = gradient;
			surface_norm_[index_i] = gradient / (gradient.norm() + TinyReal);
		}
		//=================================================================================================//
	}
//=================================================================================================//
}
//=================================================================================================//