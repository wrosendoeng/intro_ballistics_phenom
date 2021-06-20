!MODULE TO HOLD CODE PARAMETERS
!ACHILLE ARANTES BASSI, Prof
!IME 2019

module parameters
  implicit none
  integer(4),parameter:: PARM_SP = 4, PARM_DP = 8, PARM_SI = 4, PARM_LI = 8, &
    PARM_SCL = 40, PARM_LCL = 100
  
  !SP = Single Precision Real  
  !DP = Double precision real
  !SI = Integer (short)
  !LI = Integet (long)
  !WP = Working Precision
  !IntLength = Integer Length
  !SCL = Short Character Length
  !LCL = Long Character Length
    
end module parameters
  
module constants
  use parameters

  implicit none
  real(PARM_DP), parameter:: CST_PI_DP = 0.314159265358979311599796346854E+01_PARM_DP, &
    CST_E_DP = 0.271828182845904509079559829843E+01_PARM_DP
  real(PARM_SP), parameter:: CST_PI_SP = 0.314159265358979311599796346854E+01_PARM_SP, &
    CST_E_SP = 0.271828182845904509079559829843E+01_PARM_SP
end module constants