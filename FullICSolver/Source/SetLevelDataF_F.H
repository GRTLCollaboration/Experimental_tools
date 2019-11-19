#ifndef _SETLEVELDATAF_F_H_
#define _SETLEVELDATAF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGETLAPLACIANPSIF 
#define GUARDGETLAPLACIANPSIF 
// Prototype for Fortran procedure GETLAPLACIANPSIF ...
//
void FORTRAN_NAME( GETLAPLACIANPSIF ,getlaplacianpsif )(
      CHFp_FRA1(l_of_psi)
      ,CHFp_CONST_FRA1(psi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_GETLAPLACIANPSIF FORTRAN_NAME( inlineGETLAPLACIANPSIF, inlineGETLAPLACIANPSIF)
#define FORTNT_GETLAPLACIANPSIF FORTRAN_NAME( GETLAPLACIANPSIF, getlaplacianpsif)

inline void FORTRAN_NAME(inlineGETLAPLACIANPSIF, inlineGETLAPLACIANPSIF)(
      CHFp_FRA1(l_of_psi)
      ,CHFp_CONST_FRA1(psi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETLAPLACIANPSIF");
 FORTRAN_NAME( GETLAPLACIANPSIF ,getlaplacianpsif )(
      CHFt_FRA1(l_of_psi)
      ,CHFt_CONST_FRA1(psi)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETLAPLACIANPSIF 

#ifndef GUARDGETRHOGRADPHIF 
#define GUARDGETRHOGRADPHIF 
// Prototype for Fortran procedure GETRHOGRADPHIF ...
//
void FORTRAN_NAME( GETRHOGRADPHIF ,getrhogradphif )(
      CHFp_FRA1(rho_grad_phi)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_GETRHOGRADPHIF FORTRAN_NAME( inlineGETRHOGRADPHIF, inlineGETRHOGRADPHIF)
#define FORTNT_GETRHOGRADPHIF FORTRAN_NAME( GETRHOGRADPHIF, getrhogradphif)

inline void FORTRAN_NAME(inlineGETRHOGRADPHIF, inlineGETRHOGRADPHIF)(
      CHFp_FRA1(rho_grad_phi)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETRHOGRADPHIF");
 FORTRAN_NAME( GETRHOGRADPHIF ,getrhogradphif )(
      CHFt_FRA1(rho_grad_phi)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETRHOGRADPHIF 

#ifndef GUARDGETLAPLACIANUF 
#define GUARDGETLAPLACIANUF 
// Prototype for Fortran procedure GETLAPLACIANUF ...
//
void FORTRAN_NAME( GETLAPLACIANUF ,getlaplacianuf )(
      CHFp_FRA1(l_of_U)
      ,CHFp_CONST_FRA1(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_GETLAPLACIANUF FORTRAN_NAME( inlineGETLAPLACIANUF, inlineGETLAPLACIANUF)
#define FORTNT_GETLAPLACIANUF FORTRAN_NAME( GETLAPLACIANUF, getlaplacianuf)

inline void FORTRAN_NAME(inlineGETLAPLACIANUF, inlineGETLAPLACIANUF)(
      CHFp_FRA1(l_of_U)
      ,CHFp_CONST_FRA1(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETLAPLACIANUF");
 FORTRAN_NAME( GETLAPLACIANUF ,getlaplacianuf )(
      CHFt_FRA1(l_of_U)
      ,CHFt_CONST_FRA1(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETLAPLACIANUF 

#ifndef GUARDGETLAPLACIANV_IF 
#define GUARDGETLAPLACIANV_IF 
// Prototype for Fortran procedure GETLAPLACIANV_IF ...
//
void FORTRAN_NAME( GETLAPLACIANV_IF ,getlaplacianv_if )(
      CHFp_FRA(l_of_V)
      ,CHFp_CONST_FRA(V)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_GETLAPLACIANV_IF FORTRAN_NAME( inlineGETLAPLACIANV_IF, inlineGETLAPLACIANV_IF)
#define FORTNT_GETLAPLACIANV_IF FORTRAN_NAME( GETLAPLACIANV_IF, getlaplacianv_if)

inline void FORTRAN_NAME(inlineGETLAPLACIANV_IF, inlineGETLAPLACIANV_IF)(
      CHFp_FRA(l_of_V)
      ,CHFp_CONST_FRA(V)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETLAPLACIANV_IF");
 FORTRAN_NAME( GETLAPLACIANV_IF ,getlaplacianv_if )(
      CHFt_FRA(l_of_V)
      ,CHFt_CONST_FRA(V)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETLAPLACIANV_IF 

#ifndef GUARDGETGRADPHIF 
#define GUARDGETGRADPHIF 
// Prototype for Fortran procedure GETGRADPHIF ...
//
void FORTRAN_NAME( GETGRADPHIF ,getgradphif )(
      CHFp_FRA(grad_phi)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_GETGRADPHIF FORTRAN_NAME( inlineGETGRADPHIF, inlineGETGRADPHIF)
#define FORTNT_GETGRADPHIF FORTRAN_NAME( GETGRADPHIF, getgradphif)

inline void FORTRAN_NAME(inlineGETGRADPHIF, inlineGETGRADPHIF)(
      CHFp_FRA(grad_phi)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETGRADPHIF");
 FORTRAN_NAME( GETGRADPHIF ,getgradphif )(
      CHFt_FRA(grad_phi)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETGRADPHIF 

}

#endif
