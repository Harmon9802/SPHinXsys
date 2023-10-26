/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	eulerian_fluid_integration.h
 * @brief 	Here, we define the common weakly compressible eulerian classes for fluid dynamics.
 * @author	Zhentong Wang and Xiangyu Hu
 */
#ifndef EULERIAN_FLUID_INTEGRATION_H
#define EULERIAN_FLUID_INTEGRATION_H

#include "fluid_integration.hpp"
#include "riemann_solver.h"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @struct EulerianAcousticRiemannSolver
 * @brief  Acoustic RiemannSolver for Eulerian weakly-compressible flow.
 */
class EulerianAcousticRiemannSolver
{
    Fluid &fluid_i_, &fluid_j_;
    Real limiter_parameter_;

  public:
    EulerianAcousticRiemannSolver(Fluid &fluid_i, Fluid &fluid_j, Real limiter_parameter = 15.0)
        : fluid_i_(fluid_i), fluid_j_(fluid_j), limiter_parameter_(limiter_parameter){};
    FluidStarState getInterfaceState(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij);
};

template <class DataDelegationType>
class EulerianIntegration : public BaseIntegration<DataDelegationType>
{
  public:
    template <class BaseRelationType>
    explicit EulerianIntegration(BaseRelationType &base_relation);
    virtual ~EulerianIntegration(){};

  protected:
    StdLargeVec<Vecd> &mom_, &dmom_dt_;
};

template <typename... InteractionTypes>
class EulerianIntegration1stHalf;

template <class RiemannSolverType>
class EulerianIntegration1stHalf<Inner<>, RiemannSolverType>
    : public EulerianIntegration<FluidDataInner>
{
  public:
    explicit EulerianIntegration1stHalf(BaseInnerRelation &inner_relation, Real limiter_parameter = 15.0);
    template <typename BodyRelationType, typename FirstArg>
    explicit EulerianIntegration1stHalf(ConstructorArgs<BodyRelationType, FirstArg> parameters)
        : EulerianIntegration1stHalf(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~EulerianIntegration1stHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
};
using EulerianIntegration1stHalfInnerAcousticRiemann = EulerianIntegration1stHalf<Inner<>, EulerianAcousticRiemannSolver>;

template <class RiemannSolverType>
class EulerianIntegration1stHalf<ContactWall<>, RiemannSolverType>
    : public InteractionWithWall<EulerianIntegration>
{
  public:
    EulerianIntegration1stHalf(BaseContactRelation &wall_contact_relation, Real limiter_parameter = 15.0);
    template <typename BodyRelationType, typename FirstArg>
    explicit EulerianIntegration1stHalf(ConstructorArgs<BodyRelationType, FirstArg> parameters)
        : EulerianIntegration1stHalf(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~EulerianIntegration1stHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
};
using EulerianIntegration1stHalfWithWallAcousticRiemann = EulerianIntegration1stHalf<ContactWall<>, EulerianAcousticRiemannSolver>;

template <typename... InteractionTypes>
class EulerianIntegration2ndHalf;

template <class RiemannSolverType>
class EulerianIntegration2ndHalf<Inner<>, RiemannSolverType>
    : public EulerianIntegration<FluidDataInner>
{
  public:
    explicit EulerianIntegration2ndHalf(BaseInnerRelation &inner_relation, Real limiter_parameter = 15.0);
    template <typename BodyRelationType, typename FirstArg>
    explicit EulerianIntegration2ndHalf(ConstructorArgs<BodyRelationType, FirstArg> parameters)
        : EulerianIntegration2ndHalf(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~EulerianIntegration2ndHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
};
using EulerianIntegration2ndHalfAcousticRiemann = EulerianIntegration2ndHalf<Inner<>, EulerianAcousticRiemannSolver>;

/**
 * @class EulerianIntegration2ndHalfWithWall
 * @brief template density relaxation scheme with using  Riemann solver.
 */
template <class RiemannSolverType>
class EulerianIntegration2ndHalf<ContactWall<>, RiemannSolverType>
    : public InteractionWithWall<EulerianIntegration>
{
  public:
    EulerianIntegration2ndHalf(BaseContactRelation &wall_contact_relation, Real limiter_parameter = 15.0);
    template <typename BodyRelationType, typename FirstArg>
    explicit EulerianIntegration2ndHalf(ConstructorArgs<BodyRelationType, FirstArg> parameters)
        : EulerianIntegration2ndHalf(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~EulerianIntegration2ndHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
};
using EulerianIntegration2ndHalfWithWallAcousticRiemann = EulerianIntegration2ndHalf<ContactWall<>, EulerianAcousticRiemannSolver>;
} // namespace fluid_dynamics
} // namespace SPH
#endif // EULERIAN_FLUID_INTEGRATION_H
