//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Test_PETSC.h"

#include "libmesh/quadrature.h"

extern PetscErrorCode FormJacobian1(SNES, Vec, Mat, Mat, void *);
extern PetscErrorCode FormFunction1(SNES, Vec, Vec, void *);

registerMooseObject("MooseTestApp", Test_PETSC);

template <>
InputParameters
validParams<Test_PETSC>()
{
  InputParameters params = validParams<GenericConstantMaterial>();
  return params;
}

Test_PETSC::Test_PETSC(const InputParameters & parameters) : GenericConstantMaterial(parameters)
{
  int _petsc_ierr = SNESCreate(PETSC_COMM_SELF, &_petsc_snes);
  // CHKERRQ(_petsc_ierr);
  _petsc_ierr = VecCreate(PETSC_COMM_SELF, &_petsc_x);
  // CHKERRQ(_petsc_ierr);
  _petsc_ierr = VecSetSizes(_petsc_x, PETSC_DECIDE, 2);
  // CHKERRQ(_petsc_ierr);
  _petsc_ierr = VecSetFromOptions(_petsc_x);
  // CHKERRQ(_petsc_ierr);
  _petsc_ierr = VecDuplicate(_petsc_x, &_petsc_r);
  // CHKERRQ(_petsc_ierr);

  /*
     Create Jacobian matrix data structure
  */
  _petsc_ierr = MatCreate(PETSC_COMM_SELF, &_petsc_J);
  // CHKERRQ(_petsc_ierr);
  _petsc_ierr = MatSetSizes(_petsc_J, PETSC_DECIDE, PETSC_DECIDE, 2, 2);
  // CHKERRQ(_petsc_ierr);
  _petsc_ierr = MatSetFromOptions(_petsc_J);
  // CHKERRQ(_petsc_ierr);
  _petsc_ierr = MatSetUp(_petsc_J);
  // CHKERRQ(_petsc_ierr);

  /*
   Set function evaluation routine and vector.
  */
  SNESSetFunction(_petsc_snes, _petsc_r, FormFunction1, NULL);
  // CHKERRQ(_petsc_ierr);

  /*
   Set Jacobian matrix data structure and Jacobian evaluation routine
  */

  SNESSetJacobian(_petsc_snes, _petsc_J, _petsc_J, FormJacobian1, NULL);
  // CHKERRQ(_petsc_ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize nonlinear solver; set runtime options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Set linear solver defaults for this problem. By extracting the
     KSP and PC contexts from the SNES context, we can then
     directly call any KSP and PC routines to set various options.
  */
  _petsc_ierr = SNESGetKSP(_petsc_snes, &_petsc_ksp);
  // CHKERRQ(_petsc_ierr);
  _petsc_ierr = KSPGetPC(_petsc_ksp, &_petsc_pc);
  // CHKERRQ(_petsc_ierr);
  _petsc_ierr = PCSetType(_petsc_pc, PCNONE);
  // CHKERRQ(_petsc_ierr);
  _petsc_ierr = KSPSetTolerances(_petsc_ksp, 1.e-4, PETSC_DEFAULT, PETSC_DEFAULT, 20);
  // CHKERRQ(_petsc_ierr);
  /*
     Set SNES/KSP/KSP/PC runtime options, e.g.,
         -snes_view -snes_monitor -ksp_type <_petsc_ksp> -pc_type <_petsc_pc>
     These options will override those specified above as long as
     SNESSetFromOptions() is called _after_ any other customization
     routines.
  */
  _petsc_ierr = SNESSetFromOptions(_petsc_snes);
  // CHKERRQ(_petsc_ierr);
}

void
Test_PETSC::initQpStatefulProperties()
{
  computeQpProperties();
}

void
Test_PETSC::computeQpProperties()
{

  for (unsigned int i = 0; i < _num_props; i++)
    (*_properties[i])[_qp] = _prop_values[i];

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Evaluate initial guess; then solve nonlinear system
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  _petsc_ierr = VecGetArray(_petsc_x, &xx);
  // CHKERRQ(_petsc_ierr);
  xx[0] = 2.0;
  xx[1] = 3.0;
  _petsc_ierr = VecRestoreArray(_petsc_x, &xx);
  // CHKERRQ(_petsc_ierr);
  /*
     Note: The user should initialize the vector, x, with the initial guess
     for the nonlinear solver prior to calling SNESSolve().  In particular,
     to employ an initial guess of zero, the user should explicitly set
     this vector to zero by calling VecSet().
  */

  _petsc_ierr = SNESSolve(_petsc_snes, NULL, _petsc_x);
  // CHKERRQ(_petsc_ierr);
  // if (flg) {
  //   Vec f;
  //   _petsc_ierr = VecView(x, PETSC_VIEWER_STDOUT_WORLD);
  //   CHKERRQ(_petsc_ierr);
  //   _petsc_ierr = SNESGetFunction(_petsc_snes, &f, 0, 0);
  //   _petsc_ierr = PetscFinalize();
  //   return _petsc_ierr;
  //   CHKERRQ(_petsc_ierr);
  //   _petsc_ierr = VecView(r, PETSC_VIEWER_STDOUT_WORLD);
  //   CHKERRQ(_petsc_ierr);
  // }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  _petsc_ierr = SNESDestroy(&_petsc_snes);
  // CHKERRQ(_petsc_ierr);
}

/* ------------------------------------------------------------------- */
/*
   FormFunction1 - Evaluates nonlinear function, F(x).

   Input Parameters:
.  snes - the SNES context
.  x    - input vector
.  ctx  - optional user-defined context

   Output Parameter:
.  f - function vector
 */
PetscErrorCode
FormFunction1(SNES snes, Vec x, Vec f, void * ctx)
{
  PetscErrorCode ierr;
  const PetscScalar * xx;
  PetscScalar * ff;

  /*
   Get pointers to vector data.
      - For default PETSc vectors, VecGetArray() returns a pointer to
        the data array.  Otherwise, the routine is implementation dependent.
      - You MUST call VecRestoreArray() when you no longer need access to
        the array.
   */
  VecGetArrayRead(x, &xx);
  VecGetArray(f, &ff);

  /* Compute function */
  ff[0] = xx[0] * xx[0] + xx[0] * xx[1] - 3.0;
  ff[1] = xx[0] * xx[1] + xx[1] * xx[1] - 6.0;

  /* Restore vectors */
  VecRestoreArrayRead(x, &xx);
  VecRestoreArray(f, &ff);
  return 0;
}

/* ------------------------------------------------------------------- */
/*
   FormJacobian1 - Evaluates Jacobian matrix.

   Input Parameters:
.  snes - the SNES context
.  x - input vector
.  dummy - optional user-defined context (not used here)

   Output Parameters:
.  jac - Jacobian matrix
.  B - optionally different preconditioning matrix
.  flag - flag indicating matrix structure
*/
PetscErrorCode
FormJacobian1(SNES snes, Vec x, Mat jac, Mat B, void * dummy)
{
  const PetscScalar * xx;
  PetscScalar A[4];
  PetscErrorCode ierr;
  PetscInt idx[2] = {0, 1};

  /*
     Get pointer to vector data
  */
  VecGetArrayRead(x, &xx);

  /*
     Compute Jacobian entries and insert into matrix.
      - Since this is such a small problem, we set all entries for
        the matrix at once.
  */
  A[0] = 2.0 * xx[0] + xx[1];
  A[1] = xx[0];
  A[2] = xx[1];
  A[3] = xx[0] + 2.0 * xx[1];
  MatSetValues(B, 2, idx, 2, idx, A, INSERT_VALUES);

  /*
     Restore vector
  */
  VecRestoreArrayRead(x, &xx);

  /*
     Assemble matrix
  */
  MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
  if (jac != B)
  {
    MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY);
  }
  return 0;
}
