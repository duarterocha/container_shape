#define JIT_ELEMENT_SHARED_LIB
#include "jitbridge.h"

static JITFuncSpec_Table_FiniteElement_t * my_func_table;
#include "jitbridge_hang.h"


// INITIAL CONDITION 
static double ElementalInitialConditions0(const JITElementInfo_t * eleminfo, int field_index,double *_x, double *_xlagr,double t,int flag,double default_val)
{
  return default_val;
}

static double ElementalDirichletConditions(const JITElementInfo_t * eleminfo, int field_index,double *_x, double *_xlagr,double t,double default_val)
{
  if (field_index==0) // DC of field velocity_x
  {
    return 0.0; 
  }
  else if (field_index==1) // DC of field velocity_y
  {
    return 0.0; 
  }
  else if (field_index==2) // DC of field velocity_theta
  {
    return 0.0; 
  }
  else if (field_index==4) // DC of field T
  {
    return -1.0000000000000000e+00; 
  }
  return default_val;
}

static double GeometricJacobian(const JITElementInfo_t * eleminfo, const double * _x)
{
  return 6.2831853071795862e+00*_x[0];
}



JIT_API void JIT_ELEMENT_init(JITFuncSpec_Table_FiniteElement_t *functable)
{
 functable->nodal_dim=2;
 functable->lagr_dim=2;
 functable->fd_jacobian=false; 
 functable->fd_position_jacobian=false; 
 functable->debug_jacobian_epsilon = 0;
 functable->stop_on_jacobian_difference = false;
 functable->numfields_Pos=4;
 functable->fieldnames_Pos=(char **)malloc(sizeof(char*)*functable->numfields_Pos);
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_Pos,0, "coordinate_x" );
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_Pos,1, "coordinate_y" );
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_Pos,2, "lagrangian_x" );
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_Pos,3, "lagrangian_y" );
 functable->numfields_C2=3;
 functable->numfields_C2_bulk=3;
 functable->numfields_C2_basebulk=3;
 functable->fieldnames_C2=(char **)malloc(sizeof(char*)*functable->numfields_C2);
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_C2,0, "velocity_x" );
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_C2,1, "velocity_y" );
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_C2,2, "velocity_theta" );
 functable->numfields_C1=2;
 functable->numfields_C1_bulk=2;
 functable->numfields_C1_basebulk=2;
 functable->fieldnames_C1=(char **)malloc(sizeof(char*)*functable->numfields_C1);
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_C1,0, "pressure" );
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_C1,1, "T" );
 functable->dominant_space=strdup("C2");
 functable->max_dt_order=0;
 functable->moving_nodes=false;
 functable->num_res_jacs=1;
 functable->ResidualAndJacobian_NoHang=(JITFuncSpec_ResidualAndJacobian_FiniteElement *)calloc(functable->num_res_jacs,sizeof(JITFuncSpec_ResidualAndJacobian_FiniteElement));
 functable->ResidualAndJacobian=(JITFuncSpec_ResidualAndJacobian_FiniteElement *)calloc(functable->num_res_jacs,sizeof(JITFuncSpec_ResidualAndJacobian_FiniteElement));
 functable->ResidualAndJacobianSteady=(JITFuncSpec_ResidualAndJacobian_FiniteElement *)calloc(functable->num_res_jacs,sizeof(JITFuncSpec_ResidualAndJacobian_FiniteElement));
 functable->shapes_required_ResJac=(JITFuncSpec_RequiredShapes_FiniteElement_t *)calloc(functable->num_res_jacs,sizeof(JITFuncSpec_RequiredShapes_FiniteElement_t));
 functable->res_jac_names=(char**)calloc(functable->num_res_jacs,sizeof(char*));
 SET_INTERNAL_FIELD_NAME(functable->res_jac_names,0, "" );

 functable->num_Z2_flux_terms = 0;
 functable->temporal_error_scales=calloc(11,sizeof(double)); 
 functable->discontinuous_refinement_exponents=calloc(11,sizeof(double));
 functable->num_ICs=1;
 functable->IC_names=(char**)calloc(functable->num_ICs,sizeof(char*));
 functable->InitialConditionFunc=(JITFuncSpec_InitialCondition_FiniteElement*)calloc(functable->num_ICs,sizeof(JITFuncSpec_InitialCondition_FiniteElement));
 SET_INTERNAL_FIELD_NAME(functable->IC_names,0, "" );
 functable->InitialConditionFunc[0]=&ElementalInitialConditions0;
 functable->DirichletConditionFunc=&ElementalDirichletConditions;
 functable->Dirichlet_set_size=8;
 functable->Dirichlet_set=(bool *)calloc(functable->Dirichlet_set_size,sizeof(bool)); 
 functable->Dirichlet_names=(char**)calloc(functable->Dirichlet_set_size,sizeof(char*));
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,0, "" );
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,1, "coordinate_y" );
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,2, "coordinate_x" );
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,3, "velocity_x" );
 functable->Dirichlet_set[3]=true; //velocity_x
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,4, "velocity_y" );
 functable->Dirichlet_set[4]=true; //velocity_y
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,5, "velocity_theta" );
 functable->Dirichlet_set[5]=true; //velocity_theta
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,6, "pressure" );
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,7, "T" );
 functable->Dirichlet_set[7]=true; //T
 functable->integration_order=0;
 functable->GeometricJacobian=&GeometricJacobian;
 my_func_table=functable;
}
