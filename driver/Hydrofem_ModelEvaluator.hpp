// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_ModelEvaluator_Args_HPP__
#define __Hydrofem_ModelEvaluator_Args_HPP__

#include "Hydrofem_EigenFEVecMat.hpp"
#include "Hydrofem_Assembler_Base.hpp"

namespace hydrofem 
{

class ModelEvaluator
{
public:

  class InArgs
  {
  public:
    
    // set the alpha parameter
    void set_alpha(const double alpha);
    // set the beta parameter
    void set_beta(const double beta);
    // set the time
    void set_time(const double time);
    // set the time step
    void set_delta_t(const double delta_t);
    // set support for u_dot
    void set_support_u_dot(bool support_u_dot);
    // set function for u
    void set_u(const std::shared_ptr<FEVector>& u);
    // set function for u_dot
    void set_u_dot(const std::shared_ptr<FEVector>& u_dot);
    
    
    // check if time derivative is supported
    bool supports_u_dot() const;
    // check if current time is supported
    bool supports_time() const;
    // check if time step is supported
    bool supports_delta_t() const;
    
    // get the alpha param 
    double get_alpha() const;
    //
    double get_beta() const;
    //
    double get_time() const;
    //
    double get_delta_t() const;
    // get function for u
    std::shared_ptr<const FEVector> get_u() const;
    // get function for u dot
    std::shared_ptr<const FEVector> get_u_dot() const;

  private:
    
    // const
    double m_alpha;
    // const 
    double m_beta;
    // current time
    double m_time;
    // time step
    double m_delta_t;
    // true if transient
    bool m_supports_u_dot;
    // true if transient
    bool m_supports_delta_t;
    // true if transient
    bool m_supports_time;
    // global solution
    std::shared_ptr<FEVector> m_u;
    // global time derivative of solution
    std::shared_ptr<FEVector> m_u_dot;
    
  };
  
  class OutArgs
  {
  public:
    
    //
    void set_f(const std::shared_ptr<FEVector>& f);
    //
    void set_W(const std::shared_ptr<FEMatrix>& W);
    
    //
    std::shared_ptr<FEVector> get_f() const;
    //
    std::shared_ptr<FEMatrix> get_W() const;
    
  private:  
    
    // problem residual
    std::shared_ptr<FEVector> m_f;
    // problem Jacobian
    std::shared_ptr<FEMatrix> m_jac;
  };
  
  
  /** \brief Main routine to compute the out_args */  
  virtual void evalModel(const ModelEvaluator::InArgs& in_args,
                         const ModelEvaluator::OutArgs& out_args) const;
  
protected:

  std::shared_ptr<Assembler_Base> m_assembler;  
  
};

/***
 * WIP: 
 * Develop a steady ME that will wrap assembler and pass it to the nonlinear solver 
 * Develop an implicit ME that will wrap assembler and seat between assembler and time stepper
 * Develop an explicit ME that will wrap assembler and seat between assembler and time stepper
 * 
 * All the ME will take in InArgs and return OutArgs
 */

}
// end namespace hydrofem

#endif /** __Hydrofem_ModelEvaluator_Args_HPP__ */

