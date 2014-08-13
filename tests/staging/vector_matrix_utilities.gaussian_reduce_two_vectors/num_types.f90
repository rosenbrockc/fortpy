MODULE num_types
! various useful types
implicit none
private
public dp, sp, si, li

integer, parameter:: dp=selected_real_kind(15,307)
integer, parameter:: sp=selected_real_kind(6,37)
integer, parameter:: si=selected_int_kind(1) ! very short integer -10..10 range
integer, parameter:: li=selected_int_kind(18) ! Big integer -10^18..10^18 range
END MODULE
