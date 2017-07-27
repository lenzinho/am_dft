function [pc] = load_poscar(filename)
% borrowed from olle

    u=fopen(filename,'r');
    pc.header=fgetl(u);                    % read header
    pc.latpar=sscanf(fgetl(u),'%f');       % read lattice parameter
    a1=sscanf(fgetl(u),'%f %f %f');     % first basis vector
    a2=sscanf(fgetl(u),'%f %f %f');     % second basis vector
    a3=sscanf(fgetl(u),'%f %f %f');     % third basis vector
    pc.bas=[a1'; a2'; a3'];                 % construct the basis
    pc.symb=regexp(fgetl(u), '([^ \s][^\s]*)', 'match');
    pc.atomspec=sscanf(fgetl(u),repmat('%f ',1,length(pc.symb)))';
    pc.na=sum(pc.atomspec);
    coordtype=fgetl(u);
    pc.rad(1:length(pc.atomspec))=1;
    l=0;
    for i=1:length(pc.atomspec);
        for j=1:pc.atomspec(i);
            l=l+1;
            pc.r(l,:)=sscanf(fgetl(u),'%f %f %f');
            pc.species(l)=i;
            pc.mass(l)=lo_massfromsymbol(pc.symb{i});
        end
    end            
fclose(u);

a1=a1*pc.latpar;
a2=a2*pc.latpar;
a3=a3*pc.latpar;

pc.vol=abs(dot(cross(a1,a2),a3));

b1=cross(a2,a3)/pc.vol;
b2=cross(a3,a1)/pc.vol;
b3=cross(a1,a2)/pc.vol;

pc.recbas(1,:)=b1;
pc.recbas(2,:)=b2;
pc.recbas(3,:)=b3;

end

function [ mass ] = lo_massfromsymbol( s )
%LO_MASSFROMSYMBOL Summary of this function goes here
%   Detailed explanation goes here

s=strtrim(s);

mass=1;

if  strcmp(s,'Ac') 
 mass=227.028000; 
 end
 if  strcmp(s,'Ag_new') 
 mass=107.868000; 
 end
 if  strcmp(s,'Ag') 
 mass=107.868000; 
 end
 if  strcmp(s,'Ag_pv') 
 mass=107.868000; 
 end
 if  strcmp(s,'Al') 
 mass=26.982000; 
 end
 if  strcmp(s,'Am') 
 mass=243.061000; 
 end
 if  strcmp(s,'Ar') 
 mass=39.949000; 
 end
 if  strcmp(s,'As_d_GW') 
 mass=74.922000; 
 end
 if  strcmp(s,'As_d') 
 mass=74.922000; 
 end
 if  strcmp(s,'As_GW') 
 mass=74.922000; 
 end
 if  strcmp(s,'As') 
 mass=74.922000; 
 end
 if  strcmp(s,'At_d') 
 mass=209.987000; 
 end
 if  strcmp(s,'At') 
 mass=209.987000; 
 end
 if  strcmp(s,'Au_new') 
 mass=196.966000; 
 end
 if  strcmp(s,'Au') 
 mass=196.966000; 
 end
 if  strcmp(s,'Ba_sv') 
 mass=137.327000; 
 end
 if  strcmp(s,'Be') 
 mass=9.013000; 
 end
 if  strcmp(s,'Be_sv') 
 mass=9.013000; 
 end
 if  strcmp(s,'B_h') 
 mass=10.811000; 
 end
 if  strcmp(s,'Bi_d') 
 mass=208.980000; 
 end
 if  strcmp(s,'Bi') 
 mass=208.980000; 
 end
 if  strcmp(s,'Bi_pv') 
 mass=208.980000; 
 end
 if  strcmp(s,'B') 
 mass=10.811000; 
 end
 if  strcmp(s,'Br') 
 mass=79.904000; 
 end
 if  strcmp(s,'B_s') 
 mass=10.811000; 
 end
 if  strcmp(s,'Ca') 
 mass=40.078000; 
 end
 if  strcmp(s,'Ca_pv') 
 mass=40.078000; 
 end
 if  strcmp(s,'Ca_sv') 
 mass=40.078000; 
 end
 if  strcmp(s,'C_d') 
 mass=12.011000; 
 end
 if  strcmp(s,'Cd') 
 mass=112.411000; 
 end
 if  strcmp(s,'Ce_3') 
 mass=140.115000; 
 end
 if  strcmp(s,'Ce_h') 
 mass=140.115000; 
 end
 if  strcmp(s,'Ce') 
 mass=140.115000; 
 end
 if  strcmp(s,'C_GW') 
 mass=12.011000; 
 end
 if  strcmp(s,'C_h_nr') 
 mass=12.011000; 
 end
 if  strcmp(s,'C_h') 
 mass=12.011000; 
 end
 if  strcmp(s,'Cl_h') 
 mass=35.453000; 
 end
 if  strcmp(s,'Cl') 
 mass=35.453000; 
 end
 if  strcmp(s,'Co_new') 
 mass=58.933000; 
 end
 if  strcmp(s,'Co') 
 mass=58.933000; 
 end
 if  strcmp(s,'Co_sv') 
 mass=58.933000; 
 end
 if  strcmp(s,'C') 
 mass=12.011000; 
 end
 if  strcmp(s,'Cr') 
 mass=51.996000; 
 end
 if  strcmp(s,'Cr_pv_new') 
 mass=51.996000; 
 end
 if  strcmp(s,'Cr_pv') 
 mass=51.996000; 
 end
 if  strcmp(s,'Cr_sv_new') 
 mass=51.996000; 
 end
 if  strcmp(s,'Cr_sv') 
 mass=51.996000; 
 end
 if  strcmp(s,'C_s') 
 mass=12.011000; 
 end
 if  strcmp(s,'Cs_sv') 
 mass=132.900000; 
 end
 if  strcmp(s,'Cu_f') 
 mass=63.546000; 
 end
 if  strcmp(s,'Cu_new') 
 mass=63.546000; 
 end
 if  strcmp(s,'Cu') 
 mass=63.546000; 
 end
 if  strcmp(s,'Cu_pvf') 
 mass=63.546000; 
 end
 if  strcmp(s,'Cu_pv') 
 mass=63.546000; 
 end
 if  strcmp(s,'Dy_3') 
 mass=162.500000; 
 end
 if  strcmp(s,'Dy') 
 mass=162.500000; 
 end
 if  strcmp(s,'Er_2') 
 mass=167.260000; 
 end
 if  strcmp(s,'Er_3') 
 mass=167.260000; 
 end
 if  strcmp(s,'Er') 
 mass=167.260000; 
 end
 if  strcmp(s,'Eu_2') 
 mass=151.965000; 
 end
 if  strcmp(s,'Eu_3') 
 mass=151.965000; 
 end
 if  strcmp(s,'Eu_GW') 
 mass=151.965000; 
 end
 if  strcmp(s,'Eu') 
 mass=151.965000; 
 end
 if  strcmp(s,'F_d_GW') 
 mass=18.998000; 
 end
 if  strcmp(s,'Fe') 
 mass=55.847000; 
 end
 if  strcmp(s,'Fe_pv_new') 
 mass=55.847000; 
 end
 if  strcmp(s,'Fe_pv') 
 mass=55.847000; 
 end
 if  strcmp(s,'Fe_sv_h') 
 mass=55.847000; 
 end
 if  strcmp(s,'Fe_sv') 
 mass=55.847000; 
 end
 if  strcmp(s,'F_h') 
 mass=18.998000; 
 end
 if  strcmp(s,'F') 
 mass=18.998000; 
 end
 if  strcmp(s,'Fr_sv') 
 mass=223.020000; 
 end
 if  strcmp(s,'F_s') 
 mass=18.998000; 
 end
 if  strcmp(s,'Ga_d_GW') 
 mass=69.723000; 
 end
 if  strcmp(s,'Ga_d') 
 mass=69.723000; 
 end
 if  strcmp(s,'Ga_h') 
 mass=69.723000; 
 end
 if  strcmp(s,'Ga') 
 mass=69.723000; 
 end
 if  strcmp(s,'Ga_s') 
 mass=69.723000; 
 end
 if  strcmp(s,'Ga_sv_GW') 
 mass=69.723000; 
 end
 if  strcmp(s,'Gd_3') 
 mass=157.250000; 
 end
 if  strcmp(s,'Gd') 
 mass=157.250000; 
 end
 if  strcmp(s,'Ge_d3') 
 mass=72.610000; 
 end
 if  strcmp(s,'Ge_d_GW2') 
 mass=72.610000; 
 end
 if  strcmp(s,'Ge_d_GW') 
 mass=72.610000; 
 end
 if  strcmp(s,'Ge_d_GW_ref') 
 mass=72.610000; 
 end
 if  strcmp(s,'Ge_d') 
 mass=72.610000; 
 end
 if  strcmp(s,'Ge_h') 
 mass=72.610000; 
 end
 if  strcmp(s,'Ge') 
 mass=72.610000; 
 end
 if  strcmp(s,'H1.25') 
 mass=1.000000; 
 end
 if  strcmp(s,'H1.5') 
 mass=1.000000; 
 end
 if  strcmp(s,'H.5') 
 mass=1.000000; 
 end
 if  strcmp(s,'H.75') 
 mass=1.000000; 
 end
 if  strcmp(s,'He') 
 mass=4.000000; 
 end
 if  strcmp(s,'Hf') 
 mass=178.490000; 
 end
 if  strcmp(s,'Hf_pv') 
 mass=178.490000; 
 end
 if  strcmp(s,'Hf_sv_GW') 
 mass=178.490000; 
 end
 if  strcmp(s,'Hf_sv') 
 mass=178.490000; 
 end
 if  strcmp(s,'Hg') 
 mass=200.590000; 
 end
 if  strcmp(s,'H_h') 
 mass=1.000000; 
 end
 if  strcmp(s,'Ho_3') 
 mass=164.930000; 
 end
 if  strcmp(s,'Ho') 
 mass=164.930000; 
 end
 if  strcmp(s,'H') 
 mass=1.000000; 
 end
 if  strcmp(s,'In_d') 
 mass=114.820000; 
 end
 if  strcmp(s,'In') 
 mass=114.820000; 
 end
 if  strcmp(s,'I') 
 mass=126.904000; 
 end
 if  strcmp(s,'Ir') 
 mass=192.220000; 
 end
 if  strcmp(s,'K_pv') 
 mass=39.098000; 
 end
 if  strcmp(s,'Kr') 
 mass=83.800000; 
 end
 if  strcmp(s,'K_sv') 
 mass=39.098000; 
 end
 if  strcmp(s,'La') 
 mass=138.900000; 
 end
 if  strcmp(s,'La_s') 
 mass=138.900000; 
 end
 if  strcmp(s,'Li') 
 mass=7.010000; 
 end
 if  strcmp(s,'Li_sv2') 
 mass=7.010000; 
 end
 if  strcmp(s,'Li_sv') 
 mass=7.010000; 
 end
 if  strcmp(s,'Lu_3') 
 mass=174.967000; 
 end
 if  strcmp(s,'Lu') 
 mass=174.967000; 
 end
 if  strcmp(s,'Mg_new') 
 mass=24.305000; 
 end
 if  strcmp(s,'Mg') 
 mass=24.305000; 
 end
 if  strcmp(s,'Mg_pv_GW') 
 mass=24.305000; 
 end
 if  strcmp(s,'Mg_pv.old') 
 mass=24.305000; 
 end
 if  strcmp(s,'Mg_pv') 
 mass=24.305000; 
 end
 if  strcmp(s,'Mg_sv') 
 mass=24.305000; 
 end
 if  strcmp(s,'Mn') 
 mass=54.938000; 
 end
 if  strcmp(s,'Mn_pv_new') 
 mass=54.938000; 
 end
 if  strcmp(s,'Mn_pv') 
 mass=54.938000; 
 end
 if  strcmp(s,'Mn_sv') 
 mass=54.938000; 
 end
 if  strcmp(s,'Mo') 
 mass=95.940000; 
 end
 if  strcmp(s,'Mo_pv_new') 
 mass=95.940000; 
 end
 if  strcmp(s,'Mo_pv') 
 mass=95.940000; 
 end
 if  strcmp(s,'Mo_sv') 
 mass=95.940000; 
 end
 if  strcmp(s,'Na') 
 mass=22.990000; 
 end
 if  strcmp(s,'Na_pv') 
 mass=22.990000; 
 end
 if  strcmp(s,'Na_sv') 
 mass=22.990000; 
 end
 if  strcmp(s,'Nb_pv') 
 mass=92.000000; 
 end
 if  strcmp(s,'Nb_sv_new') 
 mass=92.910000; 
 end
 if  strcmp(s,'Nb_sv') 
 mass=92.910000; 
 end
 if  strcmp(s,'Nd_3') 
 mass=144.240000; 
 end
 if  strcmp(s,'Nd') 
 mass=144.240000; 
 end
 if  strcmp(s,'Ne') 
 mass=20.180000; 
 end
 if  strcmp(s,'N_h') 
 mass=14.001000; 
 end
 if  strcmp(s,'Ni_new') 
 mass=58.690000; 
 end
 if  strcmp(s,'Ni') 
 mass=58.690000; 
 end
 if  strcmp(s,'Ni_pv') 
 mass=58.690000; 
 end
 if  strcmp(s,'N') 
 mass=14.001000; 
 end
 if  strcmp(s,'Np') 
 mass=237.048000; 
 end
 if  strcmp(s,'Np_s') 
 mass=237.048000; 
 end
 if  strcmp(s,'N_s_GW') 
 mass=14.001000; 
 end
 if  strcmp(s,'N_s') 
 mass=14.001000; 
 end
 if  strcmp(s,'N_vs') 
 mass=14.001000; 
 end
 if  strcmp(s,'O_GW') 
 mass=16.000000; 
 end
 if  strcmp(s,'O_h') 
 mass=16.000000; 
 end
 if  strcmp(s,'O') 
 mass=16.000000; 
 end
 if  strcmp(s,'O_s_GW') 
 mass=16.000000; 
 end
 if  strcmp(s,'O_s') 
 mass=16.000000; 
 end
 if  strcmp(s,'Os') 
 mass=190.200000; 
 end
 if  strcmp(s,'Os_pv') 
 mass=190.200000; 
 end
 if  strcmp(s,'O_sv') 
 mass=16.000000; 
 end
 if  strcmp(s,'Pa') 
 mass=231.036000; 
 end
 if  strcmp(s,'Pa_s') 
 mass=231.036000; 
 end
 if  strcmp(s,'Pb_d') 
 mass=207.200000; 
 end
 if  strcmp(s,'Pb_d_rel2') 
 mass=207.200000; 
 end
 if  strcmp(s,'Pb_d_rel') 
 mass=207.200000; 
 end
 if  strcmp(s,'Pb') 
 mass=207.200000; 
 end
 if  strcmp(s,'Pd_new') 
 mass=106.420000; 
 end
 if  strcmp(s,'Pd') 
 mass=106.420000; 
 end
 if  strcmp(s,'Pd_pv_new') 
 mass=106.420000; 
 end
 if  strcmp(s,'Pd_pv') 
 mass=106.420000; 
 end
 if  strcmp(s,'Pd_vnew') 
 mass=106.420000; 
 end
 if  strcmp(s,'P_h') 
 mass=30.974000; 
 end
 if  strcmp(s,'Pm_3') 
 mass=146.915000; 
 end
 if  strcmp(s,'Pm') 
 mass=146.915000; 
 end
 if  strcmp(s,'Po_d') 
 mass=208.942000; 
 end
 if  strcmp(s,'Po') 
 mass=208.942000; 
 end
 if  strcmp(s,'P') 
 mass=30.974000; 
 end
 if  strcmp(s,'Pr_3') 
 mass=146.915000; 
 end
 if  strcmp(s,'Pr') 
 mass=140.907000; 
 end
 if  strcmp(s,'Pt_new') 
 mass=195.080000; 
 end
 if  strcmp(s,'Pt') 
 mass=195.080000; 
 end
 if  strcmp(s,'Pt_pv') 
 mass=195.080000; 
 end
 if  strcmp(s,'Pt_pv_ZORA') 
 mass=195.080000; 
 end
 if  strcmp(s,'Pt_ZORA') 
 mass=195.080000; 
 end
 if  strcmp(s,'Pu_h') 
 mass=244.064000; 
 end
 if  strcmp(s,'Pu') 
 mass=244.064000; 
 end
 if  strcmp(s,'Pu_s') 
 mass=244.064000; 
 end
 if  strcmp(s,'Ra_sv') 
 mass=226.025000; 
 end
 if  strcmp(s,'Rb_pv') 
 mass=85.468000; 
 end
 if  strcmp(s,'Rb_sv') 
 mass=85.468000; 
 end
 if  strcmp(s,'Re') 
 mass=186.207000; 
 end
 if  strcmp(s,'Re_pv') 
 mass=186.207000; 
 end
 if  strcmp(s,'Rh_new') 
 mass=102.906000; 
 end
 if  strcmp(s,'Rh') 
 mass=102.906000; 
 end
 if  strcmp(s,'Rh_pv_new') 
 mass=102.906000; 
 end
 if  strcmp(s,'Rh_pv') 
 mass=102.906000; 
 end
 if  strcmp(s,'Rn') 
 mass=22.017000; 
 end
 if  strcmp(s,'runelements2') 
 mass=47.880000; 
 end
 if  strcmp(s,'runelements_PBE0') 
 mass=28.085000; 
 end
 if  strcmp(s,'runelements') 
 mass=12.011000; 
 end
 if  strcmp(s,'Ru_new') 
 mass=101.070000; 
 end
 if  strcmp(s,'Ru') 
 mass=101.070000; 
 end
 if  strcmp(s,'Ru_pv_new') 
 mass=101.070000; 
 end
 if  strcmp(s,'Ru_pv') 
 mass=101.070000; 
 end
 if  strcmp(s,'Ru_sv') 
 mass=101.070000; 
 end
 if  strcmp(s,'Sb') 
 mass=121.750000; 
 end
 if  strcmp(s,'Sc') 
 mass=44.956000; 
 end
 if  strcmp(s,'Sc_sv_h') 
 mass=44.956000; 
 end
 if  strcmp(s,'Sc_sv') 
 mass=44.956000; 
 end
 if  strcmp(s,'Se_GW') 
 mass=78.960000; 
 end
 if  strcmp(s,'Se') 
 mass=78.960000; 
 end
 if  strcmp(s,'S_h') 
 mass=32.066000; 
 end
 if  strcmp(s,'Si_d_GW_nr') 
 mass=28.085000; 
 end
 if  strcmp(s,'Si_d_GW') 
 mass=28.085000; 
 end
 if  strcmp(s,'Si_h_old') 
 mass=28.085000; 
 end
 if  strcmp(s,'Si_h') 
 mass=28.085000; 
 end
 if  strcmp(s,'Si_nopc') 
 mass=28.085000; 
 end
 if  strcmp(s,'Si') 
 mass=28.085000; 
 end
 if  strcmp(s,'Si_pv_GW') 
 mass=28.085000; 
 end
 if  strcmp(s,'Si_sv_GW_nr') 
 mass=28.085000; 
 end
 if  strcmp(s,'Si_sv_GW') 
 mass=28.085000; 
 end
 if  strcmp(s,'Sm_3') 
 mass=150.360000; 
 end
 if  strcmp(s,'Sm') 
 mass=150.360000; 
 end
 if  strcmp(s,'Sn_d') 
 mass=118.710000; 
 end
 if  strcmp(s,'Sn_GW') 
 mass=118.710000; 
 end
 if  strcmp(s,'Sn') 
 mass=118.710000; 
 end
 if  strcmp(s,'S') 
 mass=32.066000; 
 end
 if  strcmp(s,'Sr_sv') 
 mass=87.620000; 
 end
 if  strcmp(s,'Ta') 
 mass=180.948000; 
 end
 if  strcmp(s,'Ta_pv') 
 mass=180.948000; 
 end
 if  strcmp(s,'Tb_3') 
 mass=158.925000; 
 end
 if  strcmp(s,'Tb') 
 mass=158.925000; 
 end
 if  strcmp(s,'Tc_new') 
 mass=98.906000; 
 end
 if  strcmp(s,'Tc') 
 mass=98.906000; 
 end
 if  strcmp(s,'Tc_pv_new') 
 mass=98.906000; 
 end
 if  strcmp(s,'Tc_pv') 
 mass=98.906000; 
 end
 if  strcmp(s,'Te') 
 mass=127.600000; 
 end
 if  strcmp(s,'Te_rel') 
 mass=127.600000; 
 end
 if  strcmp(s,'Th') 
 mass=232.039000; 
 end
 if  strcmp(s,'Th_s') 
 mass=232.030000; 
 end
 if  strcmp(s,'Ti') 
 mass=47.880000; 
 end
 if  strcmp(s,'Ti_pv') 
 mass=47.880000; 
 end
 if  strcmp(s,'Ti_sv_GW') 
 mass=47.880000; 
 end
 if  strcmp(s,'Ti_sv_h') 
 mass=47.880000; 
 end
 if  strcmp(s,'Ti_sv_new2') 
 mass=47.880000; 
 end
 if  strcmp(s,'Ti_sv_new') 
 mass=47.880000; 
 end
 if  strcmp(s,'Ti_sv') 
 mass=47.880000; 
 end
 if  strcmp(s,'Tl_d') 
 mass=204.380000; 
 end
 if  strcmp(s,'Tl') 
 mass=204.380000; 
 end
 if  strcmp(s,'Tm_3') 
 mass=168.930000; 
 end
 if  strcmp(s,'Tm') 
 mass=168.930000; 
 end
 if  strcmp(s,'U') 
 mass=238.029000; 
 end
 if  strcmp(s,'U_s') 
 mass=238.029000; 
 end
 if  strcmp(s,'V') 
 mass=50.941000; 
 end
 if  strcmp(s,'V_pv') 
 mass=50.941000; 
 end
 if  strcmp(s,'V_sv_h') 
 mass=50.941000; 
 end
 if  strcmp(s,'V_sv_new') 
 mass=50.941000; 
 end
 if  strcmp(s,'V_sv') 
 mass=50.941000; 
 end
 if  strcmp(s,'W') 
 mass=183.850000; 
 end
 if  strcmp(s,'W_pv_new') 
 mass=183.850000; 
 end
 if  strcmp(s,'W_pv') 
 mass=183.850000; 
 end
 if  strcmp(s,'Xe') 
 mass=131.294000; 
 end
 if  strcmp(s,'Yb_2_n') 
 mass=173.040000; 
 end
 if  strcmp(s,'Yb_2') 
 mass=173.040000; 
 end
 if  strcmp(s,'Yb') 
 mass=173.040000; 
 end
 if  strcmp(s,'Y_sv') 
 mass=88.906000; 
 end
 if  strcmp(s,'Zn') 
 mass=65.390000; 
 end
 if  strcmp(s,'Zn_pv') 
 mass=65.390000; 
 end
 if  strcmp(s,'Zr') 
 mass=91.224000; 
 end
 if  strcmp(s,'Zr_sv_GW') 
 mass=91.224000; 
 end
 if  strcmp(s,'Zr_sv_new') 
 mass=91.224000; 
 end
 if  strcmp(s,'Zr_sv') 
 mass=91.224000; 
 end


end

