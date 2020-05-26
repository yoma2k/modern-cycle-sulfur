module patmo_ode
contains
  subroutine fex(neq,tt,nin,dy)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    use patmo_utils
    implicit none
    integer,intent(in)::neq
    real*8,intent(in)::tt,nin(neqAll)
    real*8,intent(out)::dy(neqAll)
    real*8::d_hp(cellsNumber,speciesNumber)
    real*8::d_hm(cellsNumber,speciesNumber)
    real*8::k_hp(cellsNumber)
    real*8::k_hm(cellsNumber)
    real*8::dzz_hp(cellsNumber),dzz_hm(cellsNumber)
    real*8::kzz_hp(cellsNumber),kzz_hm(cellsNumber)
    real*8::prem(cellsNumber)
    real*8::n(cellsNumber,speciesNumber)
    real*8::dn(cellsNumber,speciesNumber)
    real*8::Tgas(cellsNumber)
    real*8::n_p(cellsNumber,speciesNumber)
    real*8::n_m(cellsNumber,speciesNumber)
    real*8::m(speciesNumber),ngas(cellsNumber)
    real*8::ngas_hp(cellsNumber),ngas_hm(cellsNumber)
    real*8::ngas_p(cellsNumber),ngas_m(cellsNumber)
    real*8::Tgas_hp(cellsNumber),Tgas_hm(cellsNumber)
    real*8::Tgas_p(cellsNumber),Tgas_m(cellsNumber)
    real*8::ngas_hpp(cellsNumber)
    real*8::ngas_hmm(cellsNumber)
    real*8::ngas_hpz(cellsNumber)
    real*8::ngas_hmz(cellsNumber)
    real*8::therm_hp(cellsNumber)
    real*8::therm_hm(cellsNumber)
    real*8::dzzh_hp(cellsNumber)
    real*8::dzzh_hm(cellsNumber)
    real*8::iTgas_hp(cellsNumber)
    real*8::iTgas_hm(cellsNumber)
    integer::i,j

    !get mass of individual species
    m(:) = getSpeciesMass()

    !roll chemistry
    do i=1,speciesNumber
      n(:,i) = nin((i-1)*cellsNumber+1:(i*cellsNumber))
    end do

    !local copy of Tgas
    Tgas(:) = nin((positionTgas-1)*cellsNumber+1:(positionTgas*cellsNumber))
    ngas(:) = nTotAll(:)

    !forward grid points
    do j=1,cellsNumber-1
      dzz_hp(j) = .5d0*(diffusionDzz(j)+diffusionDzz(j+1))
      kzz_hp(j) = .5d0*(eddyKzz(j)+eddyKzz(j+1))
      Tgas_hp(j) = .5d0*(Tgas(j)+Tgas(j+1))
      Tgas_p(j) = Tgas(j+1)
      ngas_p(j) = ngas(j+1)
      ngas_hp(j) = .5d0*(ngas(j)+ngas(j+1))
      n_p(j,:) = n(j+1,:)
    end do

    !forward grid points: boundary conditions
    dzz_hp(cellsNumber) = 0d0
    kzz_hp(cellsNumber) = 0d0
    Tgas_hp(cellsNumber) = Tgas_hp(cellsNumber-1)
    Tgas_p(cellsNumber) = Tgas_p(cellsNumber-1)
    ngas_p(cellsNumber) = ngas_p(cellsNumber-1)
    ngas_hp(cellsNumber) = ngas_hp(cellsNumber-1)
    n_p(cellsNumber,:) = n_p(cellsNumber-1,:)

    !bakcward grid points
    do j=2,cellsNumber
      dzz_hm(j) = .5d0*(diffusionDzz(j)+diffusionDzz(j-1))
      kzz_hm(j) = .5d0*(eddyKzz(j)+eddyKzz(j-1))
      Tgas_hm(j) = .5d0*(Tgas(j)+Tgas(j-1))
      Tgas_m(j) = Tgas(j-1)
      ngas_m(j) = ngas(j-1)
      ngas_hm(j) = .5d0*(ngas(j)+ngas(j-1))
      n_m(j,:) = n(j-1,:)
    end do

    !backward grid points: boundary conditions
    dzz_hm(1) = 0d0
    kzz_hm(1) = 0d0
    Tgas_hm(1) = Tgas_hm(2)
    Tgas_m(1) = Tgas_m(2)
    ngas_m(1) = ngas_m(2)
    ngas_hm(1) = ngas_hm(2)
    n_m(1,:) = n_m(2,:)

    !eqn.24 of Rimmer+Helling (2015), http://arxiv.org/abs/1510.07052
    therm_hp(:) = thermalDiffusionFactor/Tgas_hp(:)*(Tgas_p(:)-Tgas(:))
    therm_hm(:) = thermalDiffusionFactor/Tgas_hm(:)*(Tgas_m(:)-Tgas(:))
    dzzh_hp(:) = 0.5d0*dzz_hp(:)*idh2(:)
    dzzh_hm(:) = 0.5d0*dzz_hm(:)*idh2(:)
    iTgas_hp(:) = 1d0/Tgas_hp(:)
    iTgas_hm(:) = 1d0/Tgas_hm(:)
    do i=1,speciesNumber
      prem(:) = (meanMolecularMass-m(i))*gravity/kboltzmann*gridSpace(:)
      d_hp(:,i) =  dzzh_hp(:) &
          * (prem(:)*iTgas_hp(:) &
          - therm_hp(:))
      d_hm(:,i) = dzzh_hm(:) &
          * (prem(:)*iTgas_hm(:) &
          - therm_hm(:))
    end do

    k_hp(:) = (kzz_hp(:)+dzz_hp(:))*idh2(:)
    k_hm(:) = (kzz_hm(:)+dzz_hm(:))*idh2(:)

    dn(:,:) = 0d0

    dn(:,patmo_idx_COS) = &
        - krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
        - krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        + krate(:,5)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        + krate(:,6)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        - krate(:,38)*n(:,patmo_idx_COS) &
        + krate(:,46)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH) &
        + krate(:,47)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
        - krate(:,48)*n(:,patmo_idx_SH)*n(:,patmo_idx_COS) &
        - krate(:,50)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        - krate(:,51)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2)

    dn(:,patmo_idx_HO2) = &
        - krate(:,11)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
        - krate(:,22)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
        + krate(:,26)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        + krate(:,27)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
        + krate(:,56)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO) &
        + krate(:,67)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
        - krate(:,71)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
        - krate(:,72)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3)

    dn(:,patmo_idx_NO) = &
        + krate(:,18)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
        - krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO)

    dn(:,patmo_idx_N) = &
        + krate(:,33)*n(:,patmo_idx_N2) &
        + krate(:,33)*n(:,patmo_idx_N2) &
        - krate(:,78)*n(:,patmo_idx_N)*n(:,patmo_idx_N) &
        - krate(:,78)*n(:,patmo_idx_N)*n(:,patmo_idx_N)

    dn(:,patmo_idx_HSO) = &
        + krate(:,11)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
        + krate(:,14)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
        - krate(:,24)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        - krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        - krate(:,56)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO) &
        - krate(:,59)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,69)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,70)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_CO2) = &
        + krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
        - krate(:,46)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_SO3) = &
        + krate(:,22)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
        + krate(:,23)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
        + krate(:,27)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
        + krate(:,28)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,30)*n(:,patmo_idx_SO3)*n(:,patmo_idx_H2O) &
        - krate(:,42)*n(:,patmo_idx_SO3) &
        - krate(:,67)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
        - krate(:,68)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
        - krate(:,72)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3) &
        - krate(:,73)*n(:,patmo_idx_SO3) &
        + krate(:,75)*n(:,patmo_idx_H2SO4)

    dn(:,patmo_idx_H2O) = &
        + krate(:,8)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
        + krate(:,11)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
        - krate(:,30)*n(:,patmo_idx_SO3)*n(:,patmo_idx_H2O) &
        - krate(:,53)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
        - krate(:,56)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO) &
        + krate(:,75)*n(:,patmo_idx_H2SO4)

    dn(:,patmo_idx_HSO2) = &
        - krate(:,26)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        + krate(:,71)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_CO) = &
        + krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,7)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
        + krate(:,38)*n(:,patmo_idx_COS) &
        - krate(:,47)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
        - krate(:,52)*n(:,patmo_idx_CO)*n(:,patmo_idx_S)

    dn(:,patmo_idx_O2) = &
        - krate(:,5)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        + krate(:,6)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        - krate(:,13)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
        + krate(:,14)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
        + krate(:,15)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        - krate(:,16)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        - krate(:,19)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
        + krate(:,20)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
        + krate(:,23)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
        - krate(:,24)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        + krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        - krate(:,26)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        - krate(:,27)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
        - krate(:,32)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
        + krate(:,39)*n(:,patmo_idx_O3) &
        - krate(:,40)*n(:,patmo_idx_O2) &
        + krate(:,50)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        - krate(:,51)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
        + krate(:,58)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
        - krate(:,59)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        - krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        + krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        + krate(:,64)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
        - krate(:,65)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
        - krate(:,68)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
        + krate(:,69)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        - krate(:,70)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH) &
        - krate(:,70)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH) &
        + krate(:,71)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
        + krate(:,72)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3) &
        + krate(:,77)*n(:,patmo_idx_O3)

    dn(:,patmo_idx_N2) = &
        - krate(:,33)*n(:,patmo_idx_N2) &
        + krate(:,78)*n(:,patmo_idx_N)*n(:,patmo_idx_N)

    dn(:,patmo_idx_CS2) = &
        - krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        - krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,41)*n(:,patmo_idx_CS2) &
        + krate(:,48)*n(:,patmo_idx_SH)*n(:,patmo_idx_COS) &
        + krate(:,49)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_SO) = &
        + krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        + krate(:,12)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        + krate(:,13)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
        - krate(:,15)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        - krate(:,16)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        - krate(:,17)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        - krate(:,18)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
        + krate(:,19)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
        + krate(:,20)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
        + krate(:,21)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
        + krate(:,43)*n(:,patmo_idx_SO2) &
        - krate(:,45)*n(:,patmo_idx_SO) &
        - krate(:,47)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
        - krate(:,49)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        - krate(:,57)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        - krate(:,58)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
        + krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        + krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        + krate(:,62)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        + krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO) &
        - krate(:,64)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
        - krate(:,65)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
        - krate(:,66)*n(:,patmo_idx_H)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_OH) = &
        - krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
        - krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        - krate(:,8)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
        + krate(:,9)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
        + krate(:,13)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
        - krate(:,17)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        - krate(:,21)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
        + krate(:,22)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
        + krate(:,24)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        - krate(:,29)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,31)*n(:,patmo_idx_H2SO4) &
        + krate(:,31)*n(:,patmo_idx_H2SO4) &
        - krate(:,36)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        - krate(:,37)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        + krate(:,46)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH) &
        + krate(:,48)*n(:,patmo_idx_SH)*n(:,patmo_idx_COS) &
        + krate(:,53)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
        - krate(:,54)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
        - krate(:,58)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
        + krate(:,62)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        + krate(:,66)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        - krate(:,67)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
        - krate(:,69)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,74)*n(:,patmo_idx_HSO3) &
        - krate(:,76)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        - krate(:,76)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        + krate(:,81)*n(:,patmo_idx_SO2) &
        + krate(:,82)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)

    dn(:,patmo_idx_O) = &
        - krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        - krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        + krate(:,5)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        - krate(:,7)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
        - krate(:,9)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
        - krate(:,12)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        + krate(:,16)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        + krate(:,19)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
        - krate(:,28)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,32)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
        - krate(:,35)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O) &
        + krate(:,39)*n(:,patmo_idx_O3) &
        + krate(:,40)*n(:,patmo_idx_O2) &
        + krate(:,40)*n(:,patmo_idx_O2) &
        + krate(:,42)*n(:,patmo_idx_SO3) &
        + krate(:,43)*n(:,patmo_idx_SO2) &
        + krate(:,45)*n(:,patmo_idx_SO) &
        + krate(:,47)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
        + krate(:,49)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        - krate(:,50)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,52)*n(:,patmo_idx_CO)*n(:,patmo_idx_S) &
        + krate(:,54)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
        + krate(:,57)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        - krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,64)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
        + krate(:,73)*n(:,patmo_idx_SO3) &
        + krate(:,77)*n(:,patmo_idx_O3) &
        + krate(:,80)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_H2SO4) = &
        + krate(:,30)*n(:,patmo_idx_SO3)*n(:,patmo_idx_H2O) &
        - krate(:,31)*n(:,patmo_idx_H2SO4) &
        - krate(:,75)*n(:,patmo_idx_H2SO4) &
        + krate(:,76)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH)

    dn(:,patmo_idx_NO2) = &
        - krate(:,18)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
        + krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO)

    dn(:,patmo_idx_SO4) = &
        + krate(:,34)*n(:,patmo_idx_SO2) &
        - krate(:,79)*n(:,patmo_idx_SO4)

    dn(:,patmo_idx_S) = &
        + krate(:,7)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
        - krate(:,19)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
        - krate(:,20)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
        - krate(:,21)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
        + krate(:,38)*n(:,patmo_idx_COS) &
        + krate(:,41)*n(:,patmo_idx_CS2) &
        + krate(:,45)*n(:,patmo_idx_SO) &
        - krate(:,52)*n(:,patmo_idx_CO)*n(:,patmo_idx_S) &
        + krate(:,64)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
        + krate(:,65)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
        + krate(:,66)*n(:,patmo_idx_H)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_CH3SCH3) = &
        - krate(:,35)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O) &
        - krate(:,36)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        - krate(:,37)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        + krate(:,80)*n(:,patmo_idx_SO2) &
        + krate(:,81)*n(:,patmo_idx_SO2) &
        + krate(:,82)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)

    dn(:,patmo_idx_SO2) = &
        + krate(:,15)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        + krate(:,16)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        + krate(:,17)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        + krate(:,18)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
        - krate(:,22)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
        - krate(:,23)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
        + krate(:,24)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,26)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        - krate(:,28)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,29)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,31)*n(:,patmo_idx_H2SO4) &
        - krate(:,34)*n(:,patmo_idx_SO2) &
        + krate(:,35)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O) &
        + krate(:,36)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        + krate(:,37)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        + krate(:,42)*n(:,patmo_idx_SO3) &
        - krate(:,43)*n(:,patmo_idx_SO2) &
        - krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        - krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,62)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        - krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO) &
        + krate(:,67)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
        + krate(:,68)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
        - krate(:,69)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        - krate(:,71)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
        + krate(:,73)*n(:,patmo_idx_SO3) &
        + krate(:,74)*n(:,patmo_idx_HSO3) &
        - krate(:,76)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        + krate(:,79)*n(:,patmo_idx_SO4) &
        - krate(:,80)*n(:,patmo_idx_SO2) &
        - krate(:,81)*n(:,patmo_idx_SO2) &
        - krate(:,82)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)

    dn(:,patmo_idx_CH4O3S) = &
        + krate(:,37)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        - krate(:,82)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)

    dn(:,patmo_idx_HSO3) = &
        - krate(:,27)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
        + krate(:,29)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,72)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3) &
        - krate(:,74)*n(:,patmo_idx_HSO3)

    dn(:,patmo_idx_H2) = &
        + krate(:,10)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
        - krate(:,55)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_H2S) = &
        - krate(:,8)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
        - krate(:,9)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
        - krate(:,10)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
        - krate(:,11)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
        - krate(:,44)*n(:,patmo_idx_H2S) &
        + krate(:,53)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
        + krate(:,54)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
        + krate(:,55)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
        + krate(:,56)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO)

    dn(:,patmo_idx_SH) = &
        + krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
        + krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        + krate(:,8)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
        + krate(:,9)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
        + krate(:,10)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
        - krate(:,12)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        - krate(:,13)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
        - krate(:,14)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
        + krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        + krate(:,44)*n(:,patmo_idx_H2S) &
        - krate(:,46)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH) &
        - krate(:,48)*n(:,patmo_idx_SH)*n(:,patmo_idx_COS) &
        - krate(:,53)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
        - krate(:,54)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
        - krate(:,55)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
        + krate(:,57)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        + krate(:,58)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
        + krate(:,59)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        - krate(:,70)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_CS) = &
        + krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,5)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        - krate(:,6)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        - krate(:,7)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
        + krate(:,41)*n(:,patmo_idx_CS2) &
        - krate(:,49)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        + krate(:,50)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,51)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
        + krate(:,52)*n(:,patmo_idx_CO)*n(:,patmo_idx_S)

    dn(:,patmo_idx_H) = &
        - krate(:,10)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
        + krate(:,12)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        + krate(:,17)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        + krate(:,21)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
        + krate(:,44)*n(:,patmo_idx_H2S) &
        + krate(:,55)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
        - krate(:,57)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        - krate(:,62)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        - krate(:,66)*n(:,patmo_idx_H)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_O3) = &
        - krate(:,6)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        - krate(:,14)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
        - krate(:,15)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        - krate(:,20)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
        - krate(:,23)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
        - krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        + krate(:,32)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
        - krate(:,39)*n(:,patmo_idx_O3) &
        + krate(:,51)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
        + krate(:,59)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        + krate(:,65)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
        + krate(:,68)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
        + krate(:,70)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH) &
        - krate(:,77)*n(:,patmo_idx_O3)

    ngas_hpp(:) = ngas_hp(:)/ngas_p(:)
    ngas_hpz(:) = ngas_hp(:)/ngas(:)
    ngas_hmm(:) = ngas_hm(:)/ngas_m(:)
    ngas_hmz(:) = ngas_hm(:)/ngas(:)

    do i=1,chemSpeciesNumber
      dn(:,i) = dn(:,i) &
          + (k_hp(:)-d_hp(:,i)) * ngas_hpp(:) * n_p(:,i) &
          - ((k_hp(:)+d_hp(:,i)) * ngas_hpz(:) &
          + (k_hm(:)-d_hm(:,i)) * ngas_hmz(:)) * n(:,i) &
          + (k_hm(:)+d_hm(:,i)) * ngas_hmm(:) * n_m(:,i)
    end do

    !unroll chemistry
    dy(:) = 0d0
    do i=1,speciesNumber
      dy((i-1)*cellsNumber+1:(i*cellsNumber)) = dn(:,i)
    end do

  end subroutine fex
end module patmo_ode
