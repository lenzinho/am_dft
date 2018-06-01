classdef am_dft
    
    properties (Constant)
        tiny      = 1E-4; % precision of atomic coordinates
        eps       = 1E-8; % numerical precisions
        units_meV = 6.465555; % sqrt( [eV/Ang^2] * [1/amu] ) --> 6.465555 [meV]
        units_eV  = 0.006465; % sqrt( [eV/Ang^2] * [1/amu] ) --> 0.006465 [eV]
        units_THz = 9.822906; % sqrt( [eV/Ang^2] * [1/amu] ) --> 9.822906 [THz=1/ps]
        units_GHz = 9822.906; % sqrt( [eV/Ang^2] * [1/amu] ) --> 9822.906 [GHz=1/fs]
        
        amu2gram  = 1/6.02214E23;
        ev2nm     = 1239.842;
        
        usemex    = true;
        potdir    = '/Users/lenzinho/Linux/vasp.5.4/potcars/PBE.54/';
        cifdir    = '/Users/lenzinho/Developments/_materials/cif/';
    end
    
    % materials database
    
    methods (Static)
    function [uc]    = load_material(material)
        import am_dft.*
        switch material
            % toy models % NOTE: CONTINUOUS SYMMETRIES ARE NOT CHECKED FOR.
            case '1D-chain';        uc = create_cell(abc2bas([10,1],'tetra'),                  [[0;0;0]],                                                  {'H'},1);
            case '1D-dimer';        uc = create_cell(abc2bas([10,1],'tetra'),                  [[0;0;0],[0;0;0.5]],                                        {'H','He'},1);
            case '2D-square';       uc = create_cell(abc2bas([1,10],'tetra'),                  [[0;0;0],[1;1;0]/2,[1;0;0]/2,[0;1;0]/2,[0;0;1]/4],          {'H','H','H','H','H'},1); % the last atom at [0,0,1/4] breaks z-mirror symmety so that only in-plane symmetries are considered
            case '2D-BN';           uc = create_cell(abc2bas([1,1,10],'hex'),                  [[0;0;0],[2/3;1/3;0]],                                      {'B','N'},187);
            case '2D-graphene';     uc = create_cell(abc2bas([1,1,10],'hex'),                  [[2/3;1/3;0]],                                              {'C'},191);
            case '3D-cube';         uc = create_cell(abc2bas(1,'cubic'),                       [[0;0;0]],                                                  {'H'},1);
            case '3D-NaCl';         uc = create_cell(abc2bas(1,'cubic'),                       [[0;0;0], [1;1;1]/2],                                       {'Na','Cl'}, 225);
            % metals  
            case 'fcc-Co';          uc = create_cell(abc2bas(0.35441,'cubic'),                 [[0;0;0]],                                                  {'Co'},225); % ICSD 44989
            case 'hcp-Co';          uc = create_cell(abc2bas([0.25054,0.40893],'hex'),         [[1/3;2/3;1/4]],                                            {'Co'},194); % ICSD 44990
            case 'CsCl-CoFe';       uc = create_cell(abc2bas(0.28570,'cubic'),                 [[0;0;0],[1;1;1]/2],                                        {'Co','Fe'},221); % ICSD 56273
            case 'Cu';              uc = create_cell(abc2bas(0.36151,'cubic'),                 [[0;0;0]],                                                  {'Cu'},225); % ICSD 43493
            case 'fcc-Fe';          uc = create_cell(abc2bas(0.36468,'cubic'),                 [[0;0;0]],                                                  {'Fe'},225); % ICSD 44862
            case 'bcc-Fe';          uc = create_cell(abc2bas(0.29315,'cubic'),                 [[0;0;0]],                                                  {'Fe'},229); % ICSD 44863
            % salts  
            case 'NaCl';            uc = create_cell(abc2bas(0.54533,'cubic'),                 [[0;0;0], [1;1;1]/2],                                       {'Na','Cl'}, 225);
            % oxides  
            case 'Al2O3';           uc = create_cell(abc2bas([0.47617,1.29990],'hex'),         [[0;0;0.3522],[0.6936;0;0.2500]],                           {'Al','O'}, 167); % ICSD 10425
            case 'gamma-Al2O3';     uc = create_cell(abc2bas(0.79110,'cubic'),                 [[1;1;1]*0.37970,[5;5;5]/8,[0;0;0],[1;1;1]*0.15220],        {'O','Al','Al','Al'},227); % ICSD 66558 - CIF has different origin choice
            case 'eta-Al2O3';       uc = create_cell(abc2bas(0.79140,'cubic'),                 [[1;1;1]*0.37990,[5;5;5]/8,[1;0;0]*0.77390,[1;1;1]*0.19490],{'O','Al','Al','Al'},227); % ICSD 66559 - CIF has different origin choice
            case 'BiAlO3';          uc = create_cell(abc2bas([0.537546,1.33933],'hex'),        [[0;0;0],[0;0;0.2222],[0.5326;0.0099;0.9581]],              {'Bi','Al','O'},161); % ICSD 171708
            case 'Bi2Al4O9';        uc = create_cell(abc2bas([0.77134,0.81139,0.56914],'orth'),[[0.1711;0.1677;0],[0.5;0;0.2645],[0.3545;0.3399;0.5],[0;0;0.5],[0.3718;0.2056;0.2503],[0.1364;0.412;0.5],[0.1421;0.4312;0]],{'Bi','Al','Al','O','O','O','O'},55); % ICSD 88775
            case 'SrTiO3';          uc = create_cell(abc2bas(0.39010,'cubic'),                 [[0;0;0], [1;1;1]/2, [1;1;0]/2],                            {'Sr','Ti','O'}, 221); % ICSD 80873
            % case 'SrRuO3';          uc = create_cell(abc2bas([0.55729,0.78518,0.55346],'orth'),[[0;0;0],[0.5;0.25;0.99],[0.55;0.25;0.5],[0.22;0.03;0.21]], {'Ru','Sr','O','O'}, 62); % ICSD 56697
            case 'PbTiO3';          uc = create_cell(abc2bas([0.3902,0.4156],'tetra'),         [[0;0;0],[0.5;0.5;0.5377],[0.5;0.5;0.1118],[0;0.5;0.6174]], {'Pb','Ti','O','O'},99); % ICSD 61168
            case 'TbScO3';          bas = am_dft.abc2bas([5.72920,7.91700,5.46540,90,90,90]);
                                    tau = [    0.0595         0    0.4449    0.2993
                                               0.2500         0    0.2500    0.4436
                                               0.0164    0.5000    0.8798    0.3082];
                                    symb={'Tb','Sc','O','O'}; sg_code = 62;
                                    uc = am_dft.create_cell(bas,tau,symb,sg_code);
            % semiconductors
            case 'Si';              uc = create_cell(abc2bas(0.54305,'cubic'),                 [[0;0;0]],                                                  {'Si'}, 227); % ICSD 51688
            case 'GaAs';            uc = create_cell(abc2bas(0.5652,'cubic'),                  [[0;0;0], [1;1;1]/4],                                        {'Ga','As'}, 216); % ICSD 107946  
            % nitrides  
            case 'VN';              uc = create_cell(abc2bas(0.4134,'cubic'),                  [[0;0;0], [1;1;1]/2],                                       {'V','N'},  225);
            case 'ScN';             uc = create_cell(abc2bas(0.4501,'cubic'),                  [[0;0;0], [1;1;1]/2],                                       {'Sc','N'}, 225);
            case 'TiN';             uc = create_cell(abc2bas(0.4240,'cubic'),                  [[0;0;0], [1;1;1]/2],                                       {'Ti','N'}, 225);
            case 'ZrN';             uc = create_cell(abc2bas(0.4573,'cubic'),                  [[0;0;0], [1;1;1]/2],                                       {'Zr','N'}, 225);
            case 'HfN';             uc = create_cell(abc2bas(0.4524,'cubic'),                  [[0;0;0], [1;1;1]/2],                                       {'Hf','N'}, 225);
            case 'CeN';             uc = create_cell(abc2bas(0.5043,'cubic'),                  [[0;0;0], [1;1;1]/2],                                       {'Ce','N'}, 225);
            case 'CrN';             uc = create_cell(abc2bas(0.4162,'cubic'),                  [[0;0;0], [1;1;1]/2],                                       {'Cr','N'}, 225);
            otherwise; error('load_material: unknown material'); 
        end
    end 
    end
    
    % program level

    methods (Static)
        
        
        function demo_setup_calculations
            clear;clc;
            %% --------------------------------------------------------------------------------------------------
            %                          XPS Core Level Calculations
            %  --------------------------------------------------------------------------------------------------
            %% relax structures
            odir=pwd;
            dlist = {'Sn','SnO','SnO2'}; nds=numel(dlist);
            for i = 1:nds
                wdir=[odir,'/01_rlx/',dlist{i}];cd(wdir);
                    pc = am_dft.load_cell('poscar','POSCAR');
                    incar = am_dft.generate_incar('relax',{'algo','Normal'});
                    am_dft.write_incar(incar);
                    am_dft.write_poscar(pc,'POSCAR')
                    am_dft.write_kpoints(pc,[],5)
                    am_dft.write_potcar(pc);
                cd(odir);
            end
            %% XPS (Z+1 CORE LEVEL CALCULATION)
            odir=pwd;
            dlist = {'Sn','SnO','SnO2'}; nds=numel(dlist);
            for i = 1:nds
                wdir=[odir,'/01_rlx/',dlist{i}];cd(wdir);
                    pc = am_dft.load_cell('poscar','POSCAR');
                    % z+1
                    for j = 1:numel(pc.symb); if strcmp(pc.symb{j},'Sn'); pc.symb{j}='Sb'; end; end
                wdir=[odir,'/03_z+1/',dlist{i}];mkdir(wdir);cd(wdir);
                    incar = am_dft.generate_incar('scf',{'algo','Normal','ismear','-5'});
                    am_dft.write_incar(incar);
                    am_dft.write_poscar(pc,'POSCAR');
                    am_dft.write_kpoints(pc,[],5);
                    am_dft.write_potcar(pc);
                cd(odir);
            end
            %% --------------------------------------------------------------------------------------------------
            %                          What causes the dynamic instability in VN?
            %  --------------------------------------------------------------------------------------------------
            %% create pseudocubic supercell
            pc = am_dft.load_cell('poscar','POSCAR.pc');
            sc = am_dft.get_supercell(pc,'fcc');
            am_dft.write_poscar(sc,'POSCAR')
            %% relax pseudocubic structure
            pc = am_dft.load_cell('poscar','POSCAR');
            odir=pwd;wdir='0_rlx';mkdir(wdir);cd(wdir);
                incar = am_dft.generate_incar('relax',{'algo','Normal'});
                am_dft.write_incar(incar);
                am_dft.write_poscar(pc,'POSCAR')
                am_dft.write_kpoints(pc,[],5)
                am_dft.write_potcar(pc);
                am_dft.write_qsub('trio','snic2018-3-121',2,16,'','amei2@illinois.edu')
            cd(odir);
            %% scf 
            % pc = am_dft.load_cell('poscar','POSCAR');
            % odir=pwd;wdir='1_scf';mkdir(wdir);cd(wdir);
            %     incar = am_dft.generate_incar('scf',{'nbands','40','algo','Normal'}); % MBJ is not self-consistent should not be included here.
            %     am_dft.write_incar(incar);
            %     am_dft.write_poscar(pc,'POSCAR')
            %     am_dft.write_kpoints([15 15 11]);
            %     am_dft.write_potcar(pc);
            % cd(odir);
            % scf analysis
            odir=pwd;wdir='1_scf';cd(wdir); scale_energies_ = 3.5/3;
                % load all the data
                Ef         = am_dft.get_vasp('outcar:fermi');
                [uc,pc,ic] = am_dft.load_cell('poscar','POSCAR');
                [fbz,ibz]  = am_dft.get_zones(pc,[15 15 11]);
                doscar     = am_dft.load_doscar(Ef);
                dft        = am_dft.load_procar(Ef);
                ibzkt      = am_dft.load_ibzkpt(uc);
                % check that kpoints match (energies from DFT can be mapped on to the internally generated ibz points and tetrahedra)
                if sum(abs(ibz.k-ibzkt.k))>1E-8
                    error('mismatch between ibz and ibzkpt from vasp'); 
                end
                % try compling mex for faster DOS
                am_dft.compile_mex()
                % intialize figure
                set(gca,'color','w');
                am_lib.set_plot_defaults_();
                % 1) plot DOS
                subplot(2,2,1);
                dos = am_dft.get_dos(dft,ibz,doscar.E);
                plot(scale_energies_*dos.E,dos.D,'-',doscar.E,doscar.D/2,'--'); % factor of 2 for spins
                xlim([-6 6]); xlabel('E [eV]'); ylabel('g(E) [states/spin-eV]');
                % 2) compute optical joint density of states factoring in only transitions from the conduction to valence band
                subplot(2,2,2);
                jdos = am_dft.get_dos(dft,ibz,linspace(0.7,5.8,400),'jdos');
                plot(scale_energies_*jdos.E,jdos.D,'-');
                xlim([0.7 5.7]); xlabel('hv [eV]'); ylabel('jDOS(hv)');
                % 3) plot fermi surface (of conduction and valence band)
                subplot(2,2,3);
                am_dft.plot_fermi_surface(ibz,fbz,dft,[-0.3,0.3],'fermi,tetra');
                % 4) plot nesting-JDOS along high-symmetry directions
                degauss=0.1; Ep = linspace(0,6,300).';
                Aibz = plot_nesting_jdos(fbz,ibz,bzp,dft,degauss,Ep);
            cd(odir)
            %% evk
            % odir=pwd;wdir='2_evk'; copyfile('1_scf',wdir);cd(wdir);
            %     pc = am_dft.load_cell('poscar','POSCAR');
            %     incar = am_dft.generate_incar('evk',{'nbands','40','algo','Normal'});
            %     am_dft.write_incar(incar);
            %     am_dft.write_poscar(pc,'POSCAR')
            %     bzp = am_dft.get_bz_path(pc,50,'tetra');
            %     am_dft.write_kpoints(bzp);
            %     am_dft.write_potcar(pc);
            % cd(odir);
            % evk analysis

            odir=pwd;wdir='2_evk';cd(wdir);
                Ef = 6.9686; % from SCF
                [~,pc] = am_dft.load_cell('poscar','POSCAR');
                bzp= am_dft.get_bz_path(pc,50,'tetra');
                dft= am_dft.load_procar(Ef);
                % first save the image
                subplot(3,1,1:2);
                % % render the image
                % am_dft.plot_dispersion_projected(dft,bzp,pc,'orbital/atom'); ylim([-6 6]); pbaspect([2,1,1]);
                % set(gca,'XTick',[]); xlabel([]); set(gca,'YTick',[]); ylabel([]); delete(findall(gcf,'Type','line')); set(gca,'linewidth',0.1);
                % % render the image
                subplot(3,1,1:2); set(gcf,'Renderer','painters');
                am_dft.plot_dispersion_projected(dft,bzp,pc,'none'); ylim([-6 6]); pbaspect([2,1,1]); 
                delete(findall(gcf,'Type','surface')); colorbar('off');

                % then do this to save frame:
                % delete(findall(gcf,'Type','surface')); colorbar('off');
            cd(odir);
            %% opt
            % odir=pwd;wdir='3_opt'; copyfile('1_scf',wdir);cd(wdir);
            %     pc = am_dft.load_cell('poscar','POSCAR');
            %     incar = am_dft.generate_incar('opt',{'nbands','40','algo','Normal','ismear','-5'});
            %     am_dft.write_incar(incar);
            %     am_dft.write_poscar(pc,'POSCAR')
            %     am_dft.write_kpoints([15 15 11]);
            %     am_dft.write_potcar(pc);
            % cd(odir);
            % opt analysis
            am_lib.set_plot_defaults_(); scale_energies_ = 3.5/3;
            odir=pwd;wdir='3_opt';cd(wdir);
                [~,pc] = am_dft.load_cell('poscar','POSCAR');
                Ef  = am_dft.get_vasp('outcar:fermi');
                dft = am_dft.load_procar(Ef);
                [fbz,ibz] = am_dft.get_zones(pc,[15 15 11]);
                % plot e1 + i1 * e2
                set(gca,'color','w');
                ri = {'real','imag'};
                for i = 1:2
                    x = am_dft.get_vasp(['outcar:',ri{i},'_df']);
                    E = x(:,1); exx = x(:,2); exz = x(:,4);
                    subplot(2,2,i); 
                    plot(scale_energies_*E,exx,scale_energies_*E,exz); 
                    xlabel('hv [eV]'); ylabel(ri{i}); 
                    xlim([0.7 5.3]); ylim([-3 17]);
                    xticks([1:1:6]); yticks([-10:5:30]);
                    set(gca,'XMinorTick','on','YMinorTick','on')
                    legend({'exy','ez'});
                end
                % plot jdos ontop of imaginary DF
                subplot(2,2,2); 
                degauss=0.2; ex_ = x(:,1) > 0.7 & x(:,1) < 5.3; 
                jdos = am_dft.get_dos_quick(dft,ibz,E(ex_),degauss,'jdos,gauss');
                e2 = x(ex_,2); Fcv = e2.*jdos.E./jdos.D; Pcv = Fcv.*jdos.E; 
                % calculate mean
                std(Pcv(jdos.E>2))./mean(Pcv(jdos.E>2))
                hold on; plot(jdos.E,jdos.D./jdos.E.^2*100,'--'); hold off;
                legend({'exy','ez','JDOS/(hv)^2'});
                % plot experimental results (REAL)
                xyz = {'xy','z'};
                subplot(2,2,1);
                for i = 1:2
                    d=importdata(['hybrid_model/sno_',xyz{i},'.txt']);
                    E = d.data(:,1); n = d.data(:,2) + 1i*d.data(:,3);
                    hold on; h=plot(E,real(n.^2),'--'); hold off;
                end
                legend({'exy','ez','exy (exp)','ez (exp)'});
                % plot experimental results (IMAG)
                xyz = {'xy','z'};
                subplot(2,2,2);
                for i = 1:2
                    d=importdata(['hybrid_model/sno_',xyz{i},'.txt']);
                    E = d.data(:,1); n = d.data(:,2) + 1i*d.data(:,3);
                    hold on; h=plot(E,imag(n.^2),'--'); hold off;
                end
                legend({'exy','ez','JDOS','exy (exp)','ez (exp)'});
                % calculate phonon-momentum nesting function

                % plot region giving rise to the peak at 3.6 eV (color-coded based on optical matrix element)
                subplot(2,2,4);
                % load ibz optical matrix elements
                P  = am_dft.load_optical_matrix_elements();
                % average over 1,2,3
                P = sum(P(:,:,:,1,:),5)/3; 
                % map onto fbz
                P = reshape( P , [dft.nbands.^2,ibz.nks]);
                P = am_dft.ibz2fbz(fbz,ibz,P);
                P = reshape(P, [dft.nbands.^2,fbz.n]);
                % plot 
                am_dft.plot_fermi_surface(ibz,fbz,dft,3.6/scale_energies_,'jdos,tetra',P);
                am_lib.colormap('magma');
                % 
            cd(odir)
            %% eps
            odir=pwd;wdir='4_eps';copyfile('1_scf',wdir);cd(wdir);
                pc = am_dft.load_cell('poscar','POSCAR');
                incar = am_dft.generate_incar('eps',{'nbands','40','algo','Normal','ismear','-5'});
                am_dft.write_incar(incar);
                am_dft.write_poscar(pc,'POSCAR')
                am_dft.write_kpoints([15 15 11]);
                am_dft.write_potcar(pc);
            cd(odir);
            %% fermi surface (~1.5x higher density sampling, NSCF)
            % odir=pwd;wdir='6_fermi'; copyfile('1_scf',wdir); cd(wdir);
            %     pc = am_dft.load_cell('poscar','POSCAR');
            %     incar = am_dft.generate_incar('evk',{'nbands','40','algo','Normal','ismear','-5'});
            %     am_dft.write_incar(incar);
            %     am_dft.write_poscar(pc,'POSCAR')
            %     am_dft.write_kpoints([21,21,15]);
            %     am_dft.write_potcar(pc);
            % cd(odir);
            % fermi analysis
            odir=pwd;wdir='6_fermi';cd(wdir);
                % load all the data
                Ef         = am_dft.get_vasp('outcar:fermi');
                [uc,pc,ic] = am_dft.load_cell('poscar','POSCAR');
                [fbz,ibz]  = am_dft.get_zones(pc,[21,21,15]);
                dft        = am_dft.load_procar(Ef);
            %     ibzkt      = am_dft.load_ibzkpt(uc);
            %     % check that kpoints match (energies from DFT can be mapped on to the internally generated ibz points and tetrahedra)
            %     if sum(abs(ibz.k-ibzkt.k))>1E-8
            %         error('mismatch between ibz and ibzkpt from vasp'); 
            %     end
                am_dft.plot_fermi_surface(ibz,fbz,dft,[-0.3,0.3],'fermi,tetra')
                clist = am_lib.colormap_('spectral');
                colormap(clist([1,3,4],:).^(1.2));
                % plot dos
                scale_energies_ = 3.5/3;
                doscar     = am_dft.load_doscar(Ef);



            %% AIMD
            odir=pwd;wdir='7_aimd';cd(wdir);
                pc = am_dft.load_cell('poscar','POSCAR.seed');
                sc = am_dft.get_supercell(pc,eye(3)*2);
                incar = am_dft.generate_incar('aimd',{'nbands','320','algo','Normal'}); % nbands = 40*(sc.natoms/pc.natoms);
                am_dft.write_incar(incar);
                am_dft.write_poscar(sc,'POSCAR')
                am_dft.write_kpoints(1);
                am_dft.write_potcar(sc);
                am_dft.write_qsub('trio','snic2018-3-121',2,16,'','amei2@illinois.edu')
            cd(odir);
            opt analysis












            %%
            % %% --------------------------------------------------------------------------------------------------
            % %                                        New PP for Oxygen: (O_GW_new)
            % %  * turns out it is not that great. Reverting back to O_s_GW.
            % %  --------------------------------------------------------------------------------------------------
            % %% relax structure
            % pc = am_dft.load_cell('poscar','POSCAR');
            % wdir='0_rlx_pp';mkdir(wdir);cd(wdir);
            %     incar = am_dft.generate_incar('relax',{'nbands','40','algo','Normal'});
            %     am_dft.write_incar(incar);
            %     am_dft.write_poscar(pc,'POSCAR')
            %     am_dft.write_kpoints([15 15 11]);
            %     am_dft.write_potcar(pc);
            % cd(odir);
            % %% scf 
            % pc = am_dft.load_cell('poscar','POSCAR');
            % wdir='1_scf_pp';mkdir(wdir);cd(wdir);
            %     % MBJ is not self-consistent should not be included here.
            %     % also using tetrahedron integration (ismear=-5) for accurate Ef
            %     incar = am_dft.generate_incar('scf',{'nbands','40','algo','Normal','ismear','-5'});
            %     am_dft.write_incar(incar);
            %     am_dft.write_poscar(pc,'POSCAR')
            %     am_dft.write_kpoints([15 15 11]);
            %     am_dft.write_potcar(pc);
            % cd(odir);
            % %% evk (MBJ)
            % wdir='2_evk_mbj_pp'; copyfile('1_scf_pp',wdir);cd(wdir);
            %     pc = am_dft.load_cell('poscar','POSCAR');
            %     incar = am_dft.generate_incar('evk',{'nbands','40','algo','Normal','metagga','mbj'});
            %     am_dft.write_incar(incar);
            %     am_dft.write_poscar(pc,'POSCAR')
            %     bzp = am_dft.get_bz_path(pc,50,'tetra');
            %     am_dft.write_kpoints(bzp);
            %     am_dft.write_potcar(pc);
            % cd(odir);
            % %% evk (MBJ) analysis
            % wdir='2_evk_mbj_pp';cd(wdir);
            %     pc = am_dft.load_cell('poscar','POSCAR');
            %     bzp= am_dft.get_bz_path(pc,50,'tetra');
            %     Ef = am_dft.get_vasp('outcar:fermi');
            %     dft= am_dft.load_procar(Ef);
            %     am_dft.plot_dispersion_orbital_character(dft,bzp);
            %     ylim([-10 10]);
            % cd(odir); 
        end
        
        function demo_frozen_phonon_based_on_ibrion_eight()
            clear;clc;odir=pwd;
            
            % select a mode, maximum displacement, number of displacements
            modelist = [1:100]; dmax = 0.01; ndisps = 15; flag = 'A';

            % load phonon eigenvectors and unit cell
            bdir='/Users/lenzinho/Linux/calc.vasp/SnO/paper/14_frozen/DFPT';cd(bdir);
                [~,pc] = am_dft.load_cell('poscar','POSCAR');
                incar  = am_dft.load_incar('INCAR');
                d      = am_dft.get_vasp('outcar:phonon_eigendisplacements'); % obtained with iwrite = 3
                % V      = am_dft.get_vasp('outcar:phonon_eigenvectors');
                % d      = V./sqrt(pc.mass(pc.species));
                hw     = am_dft.get_vasp('outcar:phonon_energies');
            cd(odir);

            % check modelist
            modelist = modelist(modelist<=size(d,3));
            % get scale factors setting maximum displacement
            scale_ = @(d,imode) max(am_lib.normc_(d(:,:,imode)));
            % create directories
            wdir_ = @(imode,idisp) sprintf('./%.3i/%.3i',imode,idisp);
            % setup or analysis?
            switch flag
                case {'S','setup'}
                    % modify incar
                    incar.ibrion = '-1'; incar.npar  =  ''; 
                    incar.kpar   = '8';  incar.ncore = '8';
                    % set the maximum displacement to between [-0.01,0.01]; % [nm]
                    alpha = linspace(-dmax,dmax,ndisps).'/scale_(d,imode);
                    % loop
                    for i = 1:numel(modelist); imode = modelist(i);
                    for idisp = 1:numel(alpha)
                        wdir=wdir_(imode,idisp); mkdir(wdir); cd(wdir);
                            % write files
                            % POSCAR
                                sc = pc; sc.tau = sc.tau + sc.bas\(alpha(idisp)*d(:,:,imode));
                                am_dft.write_poscar(sc,'POSCAR')
                            % INCAR
                                am_dft.write_incar(incar);
                            % KPOINTS, POTCAR
                            for f = {'KPOINTS','POTCAR'}
                                copyfile([bdir,'/',f{:}],f{:});
                            end
                        cd(odir);
                    end
                    end
                case {'A','analysis'}
                    figure(1);clc;clf;set(gcf,'color','w');hold on;
                    for i = 1:numel(modelist); imode = modelist(i);
                        % grep oszicar for final energy
                        E = NaN(ndisps,1);
                        for idisp = 1:ndisps
                            wdir=wdir_(imode,idisp); if exist(wdir,'dir')~=7; continue; end; cd(wdir);
                                E(idisp) = am_dft.get_vasp('oszicar:toten');
                            cd(odir);
                        end
                        % set the maximum displacement to between [-0.01,0.01]; % [nm]
                        alpha = linspace(-dmax,dmax,ndisps).'./scale_(d,imode);
                        % plot and fit results
                        if sum(~isnan(E))>2
                            scale = max(am_lib.normc_(d(:,:,imode))); ex_ = ~isnan(E);
                            % U [eV] = 1/2 * u^2 [nm/sqrt(amu)] * amu;
                            plot(alpha(ex_),E(ex_),'o-');
                            ft_ = fit(alpha(ex_),E(ex_),'poly2');
                            plot(alpha,E,'o-',alpha,feval(ft_,alpha),'--');

                            phonon_energy = sqrt( 2 * ft_.p1 * am_dft.units_meV * 2*pi );
                            [ phonon_energy, hw(imode), phonon_energy./hw(imode) ]
                        end
                    end
            end
        end
        
        function [uc,pc,md,bvk,pp,bvt,pt] = get_phonons(opts)
            %
            % clear; clc;
            %
            % opts.continue=true;
            % opts.fposcar='POSCAR';
            % opts.fforce_position='infile.force_position';
            % opts.cutoff2=3;
            % opts.cutoff3=3; % ignored if set to zero
            % opts.dt=2; %fs
            %
            % [uc,pc,md,bvk,pp,bvt,pt] = get_phonons(opts)
            %

            import am_lib.*  am_dft.*

            % check inputs
            for f = {opts.fposcar,opts.fforce_position}
            if ~exist(f{:},'file'); error('File not found: %s',f{:}); end
            end

            % file name to save to/load from
            sname = sprintf('%s_%s_%0.2f_%0.2f.mat',...
                opts.fposcar, opts.fforce_position, opts.cutoff2, opts.cutoff3);
            if and(exist(sname,'file'),opts.continue); load(sname); else %#ok<LOAD>

                % get cells
                [uc,pc] = load_cell('poscar',opts.fposcar);
                % load md
                [md] = load_md(uc,opts.fforce_position,opts.dt);
                % get both pairs and triplets?
                if or(opts.cutoff3==0,isempty(opts.cutoff3))
                    % get pair shells
                    [bvk] = get_model('bvk', pc, opts.cutoff2, uc, md);
%                     [bvk,pp] = get_bvk(pc,uc,md,opts.cutoff2);
                    % get triplet shells
                    bvt = []; pt = [];
                else
                    % get pair shells ( do not fit for pair force constants here )
                    [bvk,pp] = get_bvk(pc,uc,md,opts.cutoff2,'-identify -model');
                    % get triplet shells ( fit for pair and triplets simultaneously here )
                    [bvt,pt,bvk] = get_bvt(pc,uc,md,opts.cutoff3,bvk,pp,'-identify -model -fit');
                end
                % save results
                save(sname,'uc','pc','md','bvk','pp','bvt','pt');
            end

            % plot results
                % plot correlation for dft vs bvk forces on atoms
                figure('color','white'); plot_bvk_vs_aimd(uc,md,bvk,pp,bvt,pt); drawnow;
        end

        function [uc,pc,dft,tb,pp] = get_tightbinding(opts)
            % clear;clc
            %
            % import am_lib.*
            
            import am_dft.*
            
            %
            % opts.continue=true;
            % opts.fposcar='POSCAR';
            % opts.feigenval='EIGENVAL.scf';
            % opts.spdf={'p','d'};
            % opts.nskips=5;
            % opts.Ef=12.8174; % grep OUTCAR for E-fermi
            % opts.cutoff2=3;
            %
            % [uc,pc,dft,tb,pp] = get_tightbinding(opts);
            %
            % NOTE: correlation coefficient for a good fit should exceed
            % 0.99. For VN with pd orbitals and a pair cutoff of 3 Ang
            % (second neighbor), R^2 = 0.998.
            %
            % defaults that should work for almost anything:
            % opts.continue=true;
            % opts.fposcar='POSCAR';
            % opts.feigenval='EIGENVAL';
            % opts.spdf=get_vasp('potcar:vrhfin:orbital');
            % opts.nskips=0;
            % opts.Ef=get_vasp('outcar:fermi');
            % opts.cutoff2=3;
            %

            import am_lib.* am_dft.*

            % check inputs
            for f = {opts.fposcar,opts.feigenval}
            if ~exist(f{:},'file'); error('File not found: %s',f{:}); end
            end

            % file name to save to/load from
            sname = sprintf('%s_%s_%02i_%0.2f_%0.2f_%s.mat',...
                opts.fposcar, opts.feigenval, opts.nskips, ...
                opts.Ef, opts.cutoff2, strjoin(opts.spdf,'+') );
            if and(exist(sname,'file'),opts.continue); load(sname); else %#ok<LOAD>

                % get cells
                [uc,pc] = load_cell('poscar',opts.fposcar);
                % load dft
                [dft]   = load_eigenval(opts.feigenval,opts.Ef);
                % construct sc
                [sc,sc.u2p,sc.p2u] = get_supercell(pc,diag([5,5,5])); sc.i2u = sc.p2u(pc.i2p); sc.u2i = pc.p2i(sc.u2p); 
                % get tb
                [tb,pp] = get_tb(pc,dft,opts.cutoff2,opts.spdf,opts.nskips);
                % save results
                save(sname,'uc','pc','dft','tb','pp');

            end

            % plot results
                % plot correlation for dft vs bvk forces on atoms
                figure('color','white'); plot_tb_vs_dft(tb,dft); drawnow;
        end

        function demo_plot_doscar()
            clear;clc;import am_dft.*
            cd('/Users/lenzinho/Linux/calc.vasp/VN/AM05.VN.ChargeDisproportionation/2_scf');
            uc = load_cell('poscar','POSCAR');
            [dos]   = load_doscar();

            % spdf={'s','p','d','f'}; k=0;
            % for l = 0:3; for m = -l:l
            % k=k+1; lm{k} = sprintf('%s%i',spdf{l+1},m);
            % end; end

            species = uc.symb(uc.species);

            % proj = [nEs,nspins,norbitals,natoms,nbands,nkpts]
            proj = permute(dos.proj,[6,1,2,3,4,5]);

            % plot
            figure(1);set(gcf,'color','w');
            [nEs,nspins,norbitals,natoms,nbands,nkpts] = size(proj); %#ok<ASGLU>
            sx = [+1,-1]; sx = sx(1:nspins);
            for i = 1:natoms
                subplot(natoms/2,2,i); 
                plot(dos.E, reshape(dos.D.*proj(:,:,:,i).*sx,nEs,[]) ); 
                title(sprintf('%i %s',i,species{i})); xlim([-10 8]);
            end

            
        end
        
        function demp_tb()
            switch 'VN'
            case 'VN'
                clear;clc; import am_dft.* am_lib.*
                cd('/Users/lenzinho/Developments/am_dft/doc/Data/VN_band_structure');
                [~,pc] = load_cell('poscar','POSCAR');
                dft = load_eigenval(get_vasp('outcar:fermi'));
                nskips = 5;
                cutoff = 0.3;
                ip = get_model('tb' , pc, cutoff, {'p','d'}, dft, nskips);
            case 'Si'
                clear;clc; import am_dft.* am_lib.*
                cd('/Users/lenzinho/Linux/calc.vasp/Si/scf');
                [~,pc] = load_cell('poscar','POSCAR');
                dft = load_eigenval(get_vasp('outcar:fermi'));
                nskips = 0;
                cutoff = 0.3;
                ip = get_model('tb' , pc, cutoff, {'sp'}, dft, nskips);
            end
        end
        
        function [varargout] = demo_toy_models(flag)
            
            switch flag
                case '1D-chain'
                    [uc,pc] = load_cell('material','1D-chain');
                    sc = get_supercell(uc,diag([1,1,5])); plot_cell(sc)
                case '1D-dimer'
                    [uc,pc] = load_cell('material','1D-dimer');
                    sc = get_supercell(uc,diag([1,1,5])); plot_cell(sc)
                case '2D-BN'
                    clear;clc;import am_dft.*
                    [uc] = load_cell('material','2D-BN');
                    sc = get_supercell(uc,diag([5,5,1])); plot_cell(sc)
                    %
                    cutoff = 1.1; spdf={'p','p'}; flag = 'toy,tb';
                    ip = get_model(flag, uc, cutoff, spdf);
                case '2D-graphene'
                    clear;clc;import am_dft.*
                    [uc,pc] = load_cell('material','2D-graphene');
                    sc = get_supercell(uc,diag([5,5,1])); plot_cell(sc)
                    %
                    cutoff = 1.1; spdf={'p'}; flag = 'toy,tb';
                    ip = get_model(flag, pc, cutoff, spdf);
                case '3D-SC'
                    % nearest neighbor tight binding on a simple cubic lattice
                    clear;clc;import am_dft.*
                    [~,pc] = load_cell('material','3D-SC');
                    cutoff = 1.1; spdf={'s'}; flag = 'toy,tb';
                    ip = get_model(flag, pc, cutoff, spdf);
                    varargout{1} = ip;
            end
            
        end
         
    end
    
    % akaiKKR
    
    methods (Static)
        
        function [uc,kkrin] = load_input()
            % c------------------------------------------------------------
            %      spc   data/vn
            % c------------------------------------------------------------
            % c   brvtyp     a        c/a   b/a   alpha   beta   gamma
            %      fcc      7.79622  ,      ,      ,      ,       ,      ,
            % c------------------------------------------------------------
            % c   edelt    ewidth    reltyp   sdftyp   magtyp   record
            %     0.001     3.5       nrl      gga91     nmag      2nd
            % c------------------------------------------------------------
            % c   outtyp    bzqlty   maxitr   pmix
            %     update      8        50    0.023
            % c------------------------------------------------------------
            % c    ntyp
            %       2
            % c------------------------------------------------------------
            % c   type    ncmp    rmt    field   mxl  anclr   conc
            %      N        2       0      0.0     1
            %                                           7      100
            %                                           0      0
            %      V        1       0      0.0     2    23     100
            % c------------------------------------------------------------
            % c   natm
            %      2
            % c------------------------------------------------------------
            % c   atmicx                        atmtyp
            %      0.0a       0.0b       0.0c     V
            %      0.5a       0.5b       0.5c     N
            % c------------------------------------------------------------
            str = load_file_(fincar);
            % 
            while true
                str{2}
                
            end
            
        end
        
        
        function [varargout]  =    get_akaikkr(flag)
            
            switch flag
                
                % DATA
                case 'data:spc_up' % spectral function E(nEs,nks), k(nEs,nks), A(nEs,nks)
                    [~,x] = system('tail -n +2 ./data/out_up.spc'); x = sscanf(x,'%f'); 
                    nEs = sum(x(:)==x(1)); x = reshape(x,3,nEs,[]); x = permute(x,[3,2,1]);
                    [varargout{1:3}] = deal(x(:,:,1).',x(:,:,2).',x(:,:,3).'); % { A, k, E }
                    
                otherwise
                    error('unknown flag');
            end
            if isempty(x); x = NaN; end
        end
        
    end
    
    % vasp
    
    methods (Static)

        % io
        
        function [bz]    = load_vasp_bz(pc,flag,varargin)
            
            feigenval='EIGENVAL';
            fprocar  ='PROCAR';
            fibzkpt  ='IBZKPT';
            fopticalmatrix = 'OPTICALMATRIX';
            
            
            if    exist(fprocar,'file')==2
                fprintf(' ... loading PROCAR'); tic;
                    [nks,nbands,natoms,norbitals,nspins,k,w,E,~,orbital,lmproj] = load_procar();
                    bz = am_bz.define(pc.bas,k,...
                        'nks',nks,'nbands',nbands,'natoms',natoms,'norbitals',norbitals,'nspins',nspins, ...
                        'w',w,'E',E,'orbital',orbital,'lmproj',lmproj);
                fprintf('(%.f secs)\n',toc);
            elseif exist(feigenval,'file')==2
                fprintf(' ... PROCAR not found. Falling back to EIGENVAL.');
                fprintf(' ... loading EIGENVAL '); tic;
                    [~,nks,nbands,k,w,E,~] = load_eigenval(feigenval);
                    bz = am_bz.define(pc.bas,k,'nks',nks,'nbands',nbands,'w',w,'E',E);
                fprintf('(%.f secs)\n',toc);
            else
                error('eigenvalues not loaded');
            end
            
            % load tetrahedra
            if    exist(fibzkpt,'file')==2
                fprintf(' ... loading IBZKPT'); tic;
                    [k,w,nks,ntets,tet,tetw] = load_ibzkpt(inv(pc.bas).'); 
                    if ~all(am_lib.eq_(nks ,bz.nks )); error('IBZKPT nks mismatch'); end
                    if ~all(am_lib.eq_(k(:),bz.k(:))); error('IBZKPT k mismatch'); end
                    if ~all(am_lib.eq_(w(:),bz.w(:))); error('IBZKPT w mismatch'); end
                    bz = bz.set('ntets',ntets,'tet',tet,'tetw',tetw);
                fprintf('(%.f secs)\n',toc);
            end
            
            % load optical matrix elements
            if exist(fopticalmatrix,'file')==2
                fprintf(' ... loading OPTICALMATRIX'); tic;
                    [nbands,~,nks,nspins,optmat] = load_optical_matrix(fopticalmatrix);
                    bz = bz.set('nbands',nbands,'nks',nks,'nspins',nspins,'optmat',optmat);
                fprintf('(%.f secs)\n',toc);
            end
            
            % shift fermi level to zero
            if exist(fopticalmatrix,'file')==2
                fprintf(' ... searhcing OUTCAR for Efermi'); tic;
                    Ef = am_dft.get_vasp('outcar:fermi');
                    bz.E = bz.E - Ef; Ef = Ef - Ef;
                fprintf('(%.f secs)\n',toc);
            else
                fprintf(' ... WARNING: OUTCAR NOT FOUND. Fermi level not shifted to zero!'); tic;
            end
                        
            cbm = min(bz.E(bz.E(:)>Ef)); [cb,ck]=am_lib.max2_(bz.E==cbm);
            vbm = max(bz.E(bz.E(:)<Ef)); [vb,vk]=am_lib.max2_(bz.E==vbm);
            Eg  = cbm-vbm;
            fprintf(' ... summary: \n');
            fprintf('     %-15s = %i\n','nks',bz.nks);
            fprintf('     %-15s = %i\n','nbands',bz.nbands);
            fprintf('     %-15s = %i\n','nspins',bz.nspins);
            fprintf('     %-15s = %i\n','natoms',bz.natoms);
            fprintf('     %-15s = %i\n','norbitals',bz.norbitals);
            fprintf('     %-15s = %+8.3f [eV] \n','fermi energy',Ef);
            fprintf('     %-15s = %+8.3f [eV] (k(:,%i) = %-5.3f,%-5.3f,%-5.3f, iband = %i)\n','CBM [eV]',cbm,ck,k(:,ck),cb);
            fprintf('     %-15s = %+8.3f [eV] (k(:,%i) = %-5.3f,%-5.3f,%-5.3f, iband = %i)\n','VBM [eV]',vbm,vk,k(:,vk),vb);
            fprintf('     %-15s = %+8.3f [eV]\n','Eg  [eV]',Eg);

            function [k,w,nks,ntets,tet,tetw] = load_ibzkpt(recbas)
                [str,~] = am_lib.load_file_('IBZKPT');
                % number of kpoints, units
                t=sscanf(str{2},'%i'); nks = t(1);
                if     contains(str{3},'R'); units = 'frac'; 
                elseif contains(str{3},'C'); units = 'cart'; 
                else;  error('units unknown'); end
                % parse points and weights
                t=sscanf(strjoin({str{[1:nks]+3}},' '),'%f'); t = reshape(t,4,[]); k = t(1:3,:); w = t(4,:); w = w./sum(w);
                % convert to fractional
                if strcmp(units,'cart'); k = recbas\k; units = 'frac'; end
                % save tetrahedra information?
                if contains(str{nks+4},'Tetra')
                   t=sscanf(str{nks+5},'%i'); ntets = t(1);
                   t=sscanf(strjoin({str{nks+5+[1:ntets]}},' '),'%i'); t = reshape(t,5,[]); tet = t(2:5,:); tetw = t(1,:);
                else
                    ntets=[];tet=[];tetw=[];
                end
            end
                
            function [nelecs,nks,nbands,k,w,E,f] = load_eigenval(feigenval)
                % load file
                str = am_lib.load_file_(feigenval);
                % parse
                buffer = sscanf(strtrim(str{6}),'%i'); % [nelecs,nks,nbands]
                [nelecs,nks,nbands]=deal(buffer(1),buffer(2),buffer(3));
                % initialize
                i=7; k=zeros(3,nks); w=zeros(1,nks); E=zeros(nbands,nks); f=zeros(nbands,nks);
                for ik = 1:nks
                    i=i+1; % skip line
                    i=i+1; buffer = sscanf(strtrim(str{i}),'%f'); % read k and weight
                    [k(1,ik),k(2,ik),k(3,ik),w(ik)]=deal(buffer(1),buffer(2),buffer(3),buffer(4));
                    % read band energies
                    for iband = 1:nbands
                        i=i+1; buffer = sscanf(strtrim(str{i}),'%f');
                        [E(iband,ik),f(iband,ik)] = deal(buffer(2),buffer(3));
                    end
                    E(:,ik) = sort(E(:,ik));
                end
            end

            function [nks,nbands,natoms,norbitals,nspins,k,w,E,f,orbital,lmproj] = load_procar()
                import am_dft.get_vasp
                nks    = am_dft.get_vasp('procar:nkpts');
                nbands = am_dft.get_vasp('procar:nbands');
                natoms = am_dft.get_vasp('procar:natoms');
                nspins = am_dft.get_vasp('procar:ispin');
                orbital= am_dft.get_vasp('procar:orbital');
                lmproj = am_dft.get_vasp('procar:lmproj');
                k = am_dft.get_vasp('procar:kpoint').';
                w = am_dft.get_vasp('procar:weight').';
                E = am_dft.get_vasp('procar:energy');
                f = am_dft.get_vasp('procar:occupation');
                norbitals = am_dft.get_vasp('procar:norbitals');
                % normalize projections by summing over atoms and characters?
                lmproj = lmproj ./ am_lib.sum_(lmproj,[2,3]);
            end

            function [nbands,nbands_cder,nks,nspins,optmat] = load_optical_matrix(fopticalmatrix) 
                % ! write to matrix elements to file
                % open(unit=123456, file="OPTICALMATRIX",action='WRITE')
                % WRITE(123456,*) WDES%NB_TOT ! , 'NBANDS'
                % WRITE(123456,*) NBANDS_CDER ! , 'NBANDS_CDER'
                % WRITE(123456,*) WDES%NKPTS  ! , 'NKS'
                % WRITE(123456,*) WDES%ISPIN  ! , 'NSPINS'
                % ! WRITE(123456,*) 'ABS(CDER_BETWEEN_STATES( NBANDS, NBANDS_CDER, NKS, ISPIN, DIR))'
                % WRITE(123456,*) ABS(CDER_BETWEEN_STATES) ! CDER_BETWEEN_STATES(WDES%NB_TOT, NBANDS_CDER, WDES%NKPTS, WDES%ISPIN, 3)
                % ! WRITE(123456,*) 'ENERGY_DER( NBANDS, NBANDS_CDER, NKS, ISPIN, DIR )'
                % WRITE(123456,*) ENERGY_DER               ! ENERGY_DER(WDES%NB_TOT, WDES%NKPTS, WDES%ISPIN, 3)
                % close(unit=123456)
                % load file
                [str,~] = am_lib.load_file_(fopticalmatrix);
                % dimensions
                buffer = str2double({str{1:4}});
                [nbands,nbands_cder,nks,nspins]=deal(buffer(1),buffer(2),buffer(3),buffer(4));
                if nbands~=nbands_cder; error('number of conduction and valence bands are not equal. I don''t know what this means.'); end
                % parse the data
                m = prod([nbands,nbands_cder,nks,nspins,3]);
                n = prod([nbands,nks,nspins,3]);
                T = sscanf(strjoin({str{5:end}},' '),'%f');
                if (m+n)~=size(T,1); error('something is wrong; either not all numbers were written or extra numbers were written that were not supposed to be'); end
                % save absolute value of optical matrix elements
                optmat = reshape( T(1:m), [nbands,nbands_cder,nks,nspins,3] );
                optmat = permute( optmat, [4,5,1,2,3] );
                
            end

        end

        function           vasp_setup(flag)
            switch flag
                case 'charge' % charge density (PARCHG)
                    odir=pwd;
                    % generate incar
                    incar = am_dft.load_incar('INCAR');
                    incar = am_dft.generate_incar(incar,...
                        {'istart','1','icharg','2','lwave','.FALSE.','lpard','.TRUE.','lsepb','.FALSE.','lsepk','.FALSE.','nbmod','-3'});
                    % these should be set manually according to each system
                    % incar.iband='28';
                    % incar.kpuse='1';
                    incar.eint='-10000'; 
                    % make directory
                    wdir=['./charge']; mkdir(wdir); cd(wdir);
                        system('ln -s ../KPOINTS KPOINTS');
                        system('ln -s ../WAVECAR WAVECAR');
                        system('ln -s ../POTCAR POTCAR');
                        system('ln -s ../POSCAR POSCAR');
                        am_dft.write_incar(incar);
                    cd(odir);
                case 'charge_noninteracting' % noninteracting charge density (CHG) : updates wavefunction but not charge density; this is probably what is desired for charge density differences
                    odir=pwd;
                    % generate incar
                    incar = am_dft.load_incar('INCAR');
                    incar = am_dft.generate_incar(incar,...
                        {'istart','0','icharg','2','lwave','.FALSE.','nelm','30','nelmdl','30','lcharg','.TRUE.'});
                    % make directory
                    wdir=['./charge_noninteracting']; mkdir(wdir); cd(wdir);
                        system('ln -s ../KPOINTS KPOINTS');
                        system('ln -s ../POTCAR POTCAR');
                        system('ln -s ../POSCAR POSCAR');
                        am_dft.write_incar(incar);
                    cd(odir);
                case 'charge_superposition' % atomic charge density (CHG) : does not update anything (wavefunction and charge density included); this is probably not what is desired for charge density differences
                    odir=pwd;
                    % generate incar
                    incar = am_dft.load_incar('INCAR');
                    incar = am_dft.generate_incar(incar,...
                        {'istart','0','icharg','12','lwave','.FALSE.','nelm','0','lcharg','.TRUE.'});
                    % make directory
                    wdir=['./charge_superposition']; mkdir(wdir); cd(wdir);
                        system('ln -s ../KPOINTS KPOINTS');
                        system('ln -s ../POTCAR POTCAR');
                        system('ln -s ../POSCAR POSCAR');
                        am_dft.write_incar(incar);
                    cd(odir);
                otherwise
                    error('unknown flag');
            end
        end
        
        function [dos]   = load_doscar(Ef)
            
            import am_lib.*
              
            [str,~] = load_file_('DOSCAR');
                % get energies, density of states, and integration DOS
                t=sscanf(str{1},'%i'); natoms = t(1);
                if exist('INCAR','file')==2
                    nEs=am_dft.get_vasp('incar:nedos'); 
                    if isnan(nEs); t=sscanf(str{6},'%f'); nEs = round(t(3)); end
                else;              t=sscanf(str{6},'%f'); nEs = round(t(3)); end
                if nargin<1; Ef = t(4); end
                nspins = round((numel(strsplit(strtrim(str{7}),' '))-1)/2);
                t=sscanf(sprintf('%s\n',str{6+[1:nEs]}),'%f'); t=reshape(t,1+2*nspins,nEs).';
                dos.E = t(:,1)-Ef; dos.D = sum(t(:,1+[1:nspins]),2); dos.iD = sum(t(:,1+nspins+[1:nspins]),2);
            % proj(nEs,nspins,norbitals,natoms) ==> proj(nspins,norbitals,natoms,nbands,nkpts,nEs)
            norbitals = (1+3+5+7);
            proj = zeros(nEs,nspins,norbitals,natoms);
            for i = [1:natoms]
                t = sscanf(sprintf('%s\n',str{6+i*(1+nEs)+[1:nEs]}),'%f');
                t = reshape(t,1+nspins*norbitals,nEs).'; % nEs x spin x orbitals
                t = reshape(t(:,2:end),nEs,nspins,norbitals);
                proj(:,:,:,i) = t;
            end
            % proj(nspins,norbitals,natoms,nbands,nkpts,nEs)
            dos.proj = permute(proj,[2,3,4,5,6,1]);
            % normalize projections
            n = sum_(dos.proj,[1,2,3,4,5]); n(eq_(n,0))=1; dos.proj=dos.proj./n; 
        end
        
        function [md]    = load_md(uc,fforces,dt)
            %
            % Loads outcar preprocessed with outcar2fp.sh (generate_script).
            %

            import am_lib.* am_dft.*

            fprintf(' ... loading displacements vs forces'); tic;

            % count number of lines in file and check that all runs completed properly
            nlines = count_lines_(fforces); if mod(nlines,uc.natoms)~=0; error('lines appear to be missing.'); end

            % open file and parse: use single precision here, solves for force constants much faster
            fid = fopen(fforces); fd = reshape(single(fscanf(fid,'%f')),6,uc.natoms,nlines/uc.natoms); fclose(fid);

            % convert [Ang] --> [uc-frac] (factor of 10 is used because uc.bas is in nm)
            fd(1:3,:,:) = matmul_(inv(uc.bas*10),fd(1:3,:,:));
            fd(4:6,:,:) = matmul_(inv(uc.bas*10),fd(4:6,:,:));

            fd_ = @(uc,force,tau,vel,dt) struct('units','frac', ...
                'bas',uc.bas,'bas2pc',uc.bas2pc,'tau2pc',uc.tau2pc,...
                'symb',{{uc.symb{:}}},'mass',uc.mass,'nspecies',uc.nspecies, ...
                'natoms',uc.natoms,'force',force,'tau',tau,'vel',vel,'species',uc.species, ...
                'dt',dt,'nsteps',size(tau,3));
            md = fd_(uc,fd(4:6,:,:),fd(1:3,:,:),cat(3,zeros(3,uc.natoms),(mod_(diff(fd(1:3,:,:),1,3)+.5)-.5)/dt),dt);

            % match to uc for saftey
            md = match_cell(md,uc);

            fprintf(' (%.f secs)\n',toc);
        end

        function [en]    = load_band_energies(nbands,fbands)
            %
            % Loads outcar preprocessed with outcar2en.sh (generate_script).
            %

            import am_lib.* am_dft.*

            % count number of lines in file and check that all runs completed properly
            nlines = count_lines_(fbands); if mod(nlines,nbands)~=0; error('lines appear to be missing.'); end

            % open file and parse
            nsteps=nlines/nbands; fid=fopen(fbands); en=reshape(fscanf(fid,'%f'),nbands,nsteps); fclose(fid);
        end

        function [incar] = load_incar(fincar)
            import am_lib.* am_dft.*
            % load incar file into memory as string
            str = load_file_(fincar);
            % generate base incar and get tokens
            incar = generate_incar('base'); token_list = fields(incar);
            for i = 1:numel(token_list)
                token = token_list{i};
                % parse
                    t = extract_token_(strtrim(str),token);
                    t = strsplit(t,'%');
                    t = strrep(t{1},'=','');
                    t = strtrim(t);
                incar.(token) = t;
            end
        end

        function incar   = generate_incar(flag,nondefault_incar_opts_)
            
            import am_dft.*
            
            if isstruct(flag)
                incar = flag;
            else
                switch flag
                    case 'base'
                        % lepsilon
                        opts_ = { ...
                        'ncore'    , '' , 'kpar'   , '' , 'npar'   , '' , 'nsim'   , '' , ... % parallelization
                        'istart'   , '' , 'icharg' , '' , 'nwrite' , '' , ... % restarting
                        'ispin'    , '' , 'gga'    , '' , 'ivdw'   , '' , ... % exchange correlations stuff
                        'encut'    , '' , 'prec'   , '' , 'algo'   , '' , 'ialgo'  , '' , 'amix'   , '' , ...
                        'ediff'    , '' , 'nelmdl' , '' , ...  % precision/convergence stuff
                        'nbands'   , '' , 'ismear' , '' , 'sigma'  , '' , 'lorbit' , '' , ...
                        'ldau'     , '' ,'ldautype', '' , 'ldaul'  , '' , 'ldauu'  , '' , 'ldauj'  , '' , ...
                        'ldauprint', '' ,'lmaxmix' , '' , ... % DFT+U stuff
                        'ibrion'   , '' , 'isif'   , '' , 'ediffg' , '' , 'potim'  , '' , 'tebeg'  , '' , ...
                        'nblock'   , '' , 'nsw'    , '' , 'smass'  , '' , 'addgrid', '' , ... % aimd stuff
                        'lreal'    , '' , 'lwave'  , '' , 'lcharg' , '' , 'lvtot'  , '' ... % save stuff
                        'lpard'    , '' , 'nbmod'  , '' , 'lsepb'  , '' , 'lsepk'  , '' , 'iband', '' , 'kpuse', '' , 'eint', '' ,  ... % charge
                        };
                        incar = generate_incar('',opts_);
                    case 'defaults'
                        opts_ = { ...
                        'istart' , '0'       , 'icharg'    , '2'       , 'gga'     , 'AM'      , 'prec'   , 'Accurate' , ...
                        'ediff'  , '1e-9'    , 'algo'      , 'VeryFast', 'encut'  , '650'      , ...
                        'ismear' , '0'       , 'sigma'     , '0.2'     , ...
                        'isif'   , '3'       , 'lorbit'    , '11'      , 'nblock'  , '1'       , 'addgrid', '.FALSE.'   , ...
                        'lreal'  , '.FALSE.' , 'lwave'     , '.FALSE.' , 'lcharg'  , '.FALSE.' , 'lvtot'  , '.FALSE.' };
                        incar = generate_incar('base',opts_);
                    case {'rlx','relax'}
                        % Setting the convergence criterion based on forces
                        % being below < 0.005 eV/Ang (based on G. Kresse's
                        % PHYSICAL REVIEW B 78, 104116 2008) works only if
                        % there is at least one ionic degree of freedom.
                        % For NaCl all ions are constrained at Wyckoff
                        % positions; therefore, ediffg = 10*ediff is more
                        % appropriate.
                        opts_ = { ...
                            'ediffg' , '1e-8'    , 'ibrion' , '2'       , 'potim'  , '0.5'    , ...
                            'nsw'    , '100'     , 'smass'  , '0'       , 'addgrid', '.TRUE.' , ...
                        };
                        incar = generate_incar('defaults',opts_);
                    case 'scf'
                        opts_ = { ...
                            'ibrion' , '-1'       , 'nwrite' , '3'       , 'ismear' , '-5', ...
                            'nsw'    , '0'        , 'lorbit' , '11'      , 'nedos'  ,  '10001' ...
                            'addgrid', '.FALSE.'  , ... % takes too long with addgrid
                            'lwave'  , '.TRUE.'   , 'lcharg' , '.TRUE.'  , 'lvtot'  , '.TRUE.'  , ...
                            };
                        incar = generate_incar('defaults',opts_);
                    case {'nscf','evk'}
                        opts_ = {'icharg','11','nwrite','3','lwave','.FALSE.','lcharg','.FALSE.','lvtot','.FALSE.'};
                        incar = generate_incar('scf',opts_);
                    case 'opt'
                        opts_ = {'loptics','.TRUE.','nedos','10001'}; 
                        incar = generate_incar('nscf',opts_);
                    case 'eps' % density functional perturbation theory (requires the calculation to be self-consistent)
                        % nsw = 1 is required otherwise ibrion = 8 option is not accepted.
                        % icharg = 2 is required for ibrion = 8 other nonorthogonality error.
                        opts_ = {'icharg','2','ibrion','8','nsw','1','lepsilon','.TRUE.','isym','0','#ncore','8'}; 
                        incar = generate_incar('nscf',opts_);
                    case 'elastic' % finite difference calculation of phonons at gamma and elastic constants;  (requires the calculation to be self-consistent)
                        % nsw = 1 is required otherwise ibrion = 6 option is not accepted.
                        % icharg = 2 is required for ibrion = 8 other nonorthogonality error.
                        % isym = 0 makes the calculation much faster because it allows ncore to be set (despite needing to consider more degrees of freedom)
                        opts_ = {'icharg','2','ibrion','6','nsw','1','isif','3','potim','0.01','isym','0','#ncore','8'}; 
                        incar = generate_incar('nscf',opts_);
                    case 'charge'
                        opts_ = {'istart','1','icharg','2','lwave','.FALSE.', ...
                            'lpard','.TRUE.','nbmod','-3','lsepb','.TRUE.','lsepk','.TRUE.', ...
                            'iband','1','kpuse','1','eint','-10000'}; % these options are material specific!
                        incar = generate_incar('nscf',opts_);
                    case 'charge_noninteracting' % noninteracting charges as reference for differential charge density
                        opts_ = {'istart','0','icharg','2','lwave','.FALSE.','nelm','30','nelmdl','30','lcharg','.TRUE.'};
                        incar = generate_incar('nscf',opts_);
                    case 'charge_superposition' % superposition of atomic charges as reference for differential charge density
                        opts_ = {'istart','0','icharg','12','lwave','.FALSE.','nelm','0','lcharg','.TRUE.'};
                        incar = generate_incar('nscf',opts_);
                    case 'aimd' 
                        % similar to rlx but with ibrion=0, a lower energy convergence, no convergence on ionic motion, nose thermostat, and more steps 
                        opts_ = {'ibrion','0','ediff','1e-6','ediffg','','smass','0','nsw', '1000','tebeg','300',};
                        incar = generate_incar('rlx',opts_);
                    case 'icore1'
                        opts_ = {'icorelevel','1'};
                        incar = generate_incar('scf',opts_);
                    case 'icore2' % this option is material specific!
                        opts_ ={'icorelevel','2', ...
                                'clnt' , '1' , ... % species which is excited
                                'cln'  , '3' , ... % main quantum number of excited core electron 
                                'cll'  , '2' , ... % l quantum number of excited core electron
                                'clz'  , '1'       % (0.5) slaters transition state , (1.0) whole electron excited 
                        };
                        incar = generate_incar('scf',opts_);
                    case {'lobster'}
                        opts_ = {'isym','0'}; % lobster requires all kpoints without symmeterization
                        incar = generate_incar('scf',opts_);
                end
            end
            
            if nargin > 1
                n = numel(nondefault_incar_opts_);
                if n > 1; for i = 1:2:n
                    incar.(nondefault_incar_opts_{i})=nondefault_incar_opts_{i+1};
                end; end
            end
            
        end

        function matrix_elements = load_wavder()
            % terrible.. most likely this will not work in every case
            % use hex fiend to look at the binary file. compared the garbage to the file i wrote explicitly.
            fname = 'WAVEDER';

            %       IF (IU0>=0) THEN
            %          OPEN(UNIT=IU,FILE=DIR_APP(1:DIR_LEN)//'WAVEDER', &
            %               FORM='UNFORMATTED',STATUS='UNKNOWN')
            %          
            %          WRITE(IU) WDES%NB_TOT, NBANDS_CDER, WDES%NKPTS, WDES%ISPIN
            %          WRITE(IU) NODES_IN_DIELECTRIC_FUNCTION
            %          WRITE(IU) WPLASMON
            %          WRITE(IU) CDER_BETWEEN_STATES
            %          CLOSE(IU)
            %       ENDIF

            fid=fopen(fname, 'rb'); % Open the file.
            fread(fid, 1, 'int32'); % record headers
            nbands = fread(fid, 1, 'int32');
            nbands_cder = fread(fid, 1, 'int32');
            nks = fread(fid, 1, 'int32');
            nspins = fread(fid, 1, 'int32');
            nnodes = fread(fid, 1, 'int32');
            fread(fid, 25, 'int32'); % garbage
            matrix_elements = fread(fid, prod([2,nbands,nbands_cder,nks,nspins,3]), 'single');
            matrix_elements = reshape(matrix_elements,[2,nbands,nbands_cder,nks,nspins,3]);
            matrix_elements = matrix_elements(1,:,:,:,:,:)+1i*matrix_elements(1,:,:,:,:,:);
            fclose(fid);
        end
        
        
        % for automating vasp simulations

        function           write_poscar(uc,fposcar)

            import am_lib.* am_dft.*

            % write_poscar(uc,fposcar)
            if nargin < 2; fposcar='POSCAR'; end

            % check that the basis is right-handed
            if det(uc.bas)<0; error('vasp requires basis vectors with positive determinants'); end

            % convert from [nm] to [Ang]
            uc.bas = uc.bas*10;
            
            n = size(uc.tau,3);
            for i = 1:n
                % file name
                fname = sprintf('%s',fposcar);
                if n ~= 1; fname = sprintf('%s%s_%06i',fname,fposcar,i); end
                % open file
                fid=fopen(fname,'w');
                    % header
                    fprintf(fid,'%s',uc.get_formula());
                    if n ~= 1; fprintf(fid,' %i of %i',i,n); end
                    fprintf(fid,'\n');
                    % print body
                    fprintf(fid,'%12.8f \n',1.0);  % latpar
                    fprintf(fid,'%12.8f %12.8f %12.8f \n',uc.bas(:,1));
                    fprintf(fid,'%12.8f %12.8f %12.8f \n',uc.bas(:,2));
                    fprintf(fid,'%12.8f %12.8f %12.8f \n',uc.bas(:,3));
                    fprintf(fid,' %s ',uc.type(:).symb); fprintf(fid,'\n');
                    fprintf(fid,' %i ',uc.nspecies); fprintf(fid,'\n');
                    fprintf(fid,'Direct \n');
                    fprintf(fid,'%12.8f %12.8f %12.8f \n',uc.tau(:,:,i));
                fclose(fid);
            end
        end

        function           write_potcar(uc,potdir)
            import am_dft.*

            % default
            if nargin<2;potdir=am_dft.potdir; end
            % database
            potcar_database = {...
                '  ' ,'  ' ,'Li_sv_GW' ,'Be_sv_GW' ,'  ' ,'  ' ,'N_s_GW' ,'O_s_GW' , ... %  h     he    li    be    b     c     n     o
                '  ' ,'  ' ,'Na_sv_GW' ,'Mg_sv_GW' ,'Al_sv_GW' ,'Si_sv_GW' ,'  ' ,'  ' , ... %  f     ne    na    mg    al    si    p     s
                '  ' ,'  ' ,'K_sv' ,'Ca_sv' ,'Sc_sv_GW' ,'Ti_pv' ,'V_sv_GW' ,'  ' , ... %  cl    ar    k     ca    sc    ti    v     cr
                'Mn_GW' ,'Fe_GW' ,'  ' ,'Ni_sv_GW' ,'  ' ,'  ' ,'Ga_sv_GW' ,'  ' , ... %  mn    fe    co    ni    cu    zn    ga    ge
                '  ' ,'  ' ,'  ' ,'  ' ,'Rb_sv' ,'Sr_sv' ,'Y_sv' ,'  ' , ... %  as    se    br    kr    rb    sr    y     zr
                '  ' ,'  ' ,'  ' ,'  ' ,'Rh_GW' ,'  ' ,'  ' ,'  ' , ... %  nb    mo    tc    ru    rh    pd    ag    cd
                'In_sv_GW' ,'Sn_sv_GW' ,'Sb_sv_GW' ,'Te_sv_GW' ,'  ' ,'  ' ,'Cs_sv_GW' ,'Ba_sv_GW' , ... %  in    sn    sb    te    i     xe    cs    ba
                'La_GW' ,'Ce_3' ,'Pr_3' ,'Nd_3' ,'Pm_3' ,'Sm_3' ,'Eu_3' ,'Gd_3' , ... %  la    ce    pr    nd    pm    sm    eu    gd
                'Tb_3' ,'Dy_3' ,'Ho_3' ,'Er_3' ,'Tm_3' ,'Yb_3' ,'Lu' ,'  ' , ... %  tb    dy    ho    er    tm    yb    lu    hf
                '  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' , ... %  ta    w     re    os    ir    pt    au    hg
                '  ' ,'  ' ,'Bi_GW' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' , ... %  tl    pb    bi    po    at    rn    fr    ra
                '  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' , ... %  ac    th    pa    u     np    pu    am    cm
                '  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' , ... %  bk    cf    es    fm    md    no    lr    rf
                '  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' ,'  ' };    %  db    sg    bh    hs    mt    ds    rg    uub
            Z = [uc.type(:).Z]; nZs = numel(Z);
            for i = 1:nZs
                if isempty(potcar_database{Z(i)}); error('potcar not defined'); end
                if i == 1
                    cmd = sprintf('cat %s/%s/POTCAR > POTCAR' ,potdir,potcar_database{Z(i)});
                else
                    cmd = sprintf('cat %s/%s/POTCAR >> POTCAR',potdir,potcar_database{Z(i)});
                end
                system(cmd);
            end
        end

        function           write_incar(incar,fincar)
            if nargin < 2; fincar='INCAR'; end
            token_list = fields(incar); fclose('all');
            fid = fopen(fincar,'w');
                for i = 1:numel(token_list)
                    token = token_list{i};
                    if ~isempty(incar.(token))
                        fprintf(fid,'%-10s ',token);
                        fprintf(fid,'= %-10s',incar.(token));
                        fprintf(fid,' %% %-10s',get_comments_(token));
                        fprintf(fid,'\n');
                    end
                end
            fclose(fid);

            function c = get_comments_(token)
                switch token
                % parallelization
                case {'ncore','#ncore'}; c='orbital parallelization SQRT( number of cores )';
                case 'npar';   		c='band parallelization';
                case 'kpar';        c='kpoint parallelization';
                case 'nsim';   		c='number of bands simultaneously optimized by RMMS-DIIS algorithm';
                % starting
                case 'istart'; 		c='(0) new  (1) cont (2) samecut';
                case 'icharg'; 		c='(1) file (2) atom (3) const <scf | nscf> (11)-chgcar';
                case 'nwrite';  	c='(2) default (3) verbose (4) debugging';
                % electronic
                case 'gga';    		c='(AM05) exchannge correlation function; not usable with ivdw';
                case 'ivdw';   		c='(0) no vdW correction (2) Tkatchenko-Scheffler; not usable with density functional perturbation theory';
                case 'nelm'; 		c='number of electronic iterations';
                case 'nelmdl'; 		c='non-selfconsistent iterations: (<0) only once, (>0) at every ionic step';
                case 'ispin';  		c='(1) non-spin-polarized (2) spin-polarized';
                case 'prec';   		c='Accurate, Normal';
                case 'algo';   		c='(Veryfast) RMM-DIIS (Fast) Davidson/RMM-DIIS (Normal) Davidson';
                case 'ialgo';  		c='(38) Davidson (48) RMM-DIIS algorithm';
                case 'nbands'; 		c='number of bands';
                case 'encut';  		c='cutoff';
                case 'ismear'; 		c='(0) gauss (1) mp (-5) tetrah';
                case 'sigma';  		c='maximize with entropy below 0.1 meV / atom';
                case 'ediff';  		c='toten convergence';
                case 'isym';        c='(0) symmetry off (1) symmetry on (2) efficient symmetry for PAW';
                % magnetic
                case 'nupdown';     c='(-1) full spin relaxation (0) paramagnetic (other) value to constrain net occupancy (up - down)';
                case 'magmom';      c='magnetic moment on each atom';
                % dft+U
                case 'ldau'; 		c='Activates DFT+U calculation';
                case 'ldautype'; 	c='(1) Liechtenstein (2) Dudarev Ueff = U - J (3) adds Ueff on the atomic sites';
                case 'ldaul'; 	    c='one entry per species: (>0) quantum number to which U is applied (-1) U not applied to this species';
                case 'ldauu'; 	    c='(for LDAUTYPE = 1,2) Hubbard U per species  (for LDAUTYPE = 3) Ueff acting on spin up';
                case 'ldauj'; 	    c='(for LDAUTYPE = 1,2) Exchange J per species (for LDAUTYPE = 3) Ueff acting on spin down';
                case 'ldauprint';   c='(0) silent (1) write occupancy matrix OUTCAR (2) occupancy matrix OUTCAR and potential matrix STDOUT';
                case 'lmaxmix'; 	c='(4) d-electrons, (6) f-electrons';
                % ionic
                case 'ediffg';      c='(<0) relaxation will stop if all forces are less than |EDIFFG| (>0) relaxation stops if energy change is less than EDIFFG';
                case 'tebeg';  		c='temperature of the MD';
                case 'ibrion'; 		c='(-1) no update (0) MD (1) RMM-DIIS quasi-Newton (2) conjugate-gradient (6) finite difference (9) pertubation theory';
                case 'potim';  		c='Timestep in femtoseconds';
                case 'isif';   		c='Relax (2) ions (3) ions,volume,shape (4) ions,shape';
                case 'nblock'; 		c='write after every iteration';
                case 'nsw';    		c='number of ionic steps';
                case 'smass';  		c='(-3) NVE (=>0) Nose thermostat';
                case 'addgrid';     c='(.TRUE.) adds support grid for augmentation charges reducing noise on force (cannot be used with ivdw)';
                % icore
                case 'icorelevel';  c='(1) initial state approximation, (2) electron is removed from the core and placed into the valence' ;
                case 'clnt';        c='(icorelevel=2): species which is excited';
                case 'cln';         c='(icorelevel=2): main quantum number of excited core electron';
                case 'cll';         c='(icorelevel=2): l quantum number of excited core electron';
                case 'clz';         c='(icorelevel=2): (0.5) slaters transition state , (1.0) whole electron excited';
                % charge
                case {'lpard','#lpard'}; c='(.TRUE.) calculate charge densities';
                case {'lsepb','#lsepb'}; c='(.TRUE.) seperate charge by bands';
                case {'lsepk','#lsepk'}; c='(.TRUE.) seperate charge by kpoints';
                case {'iband','#iband'}; c='indices of bands to use in charge density calculation';
                case {'kpuse','#kpuse'}; c='indices of kpoints to use in charge density calculation';
                case 'nbmod';            c='(-1) default (-3) relative to fermi energy';
                case {'eint','#eint'};   c='energy range for charge density calculation';
                % oribtal info
                case 'lorbit'; 		c='(11) lm-proj DOS and PROCAR (0) DOS';
                % projection scheme
                case 'lreal';  		c='(.FALSE.) for <20 atoms (Auto) for >20 atoms';
                % save
                case 'lwave';  		c='write wavefunctions';
                case 'lcharg'; 		c='write charge density';
                case 'lvtot';  		c='write local charge potential';
                otherwise;     		c='';
                end
            end
        end

        function           write_kpoints(bz,fkpoints,kdensity)
            if nargin < 2 || isempty(fkpoints); fkpoints='KPOINTS'; end
            if     isnumeric(bz)
                % check format
                if ~isvector(bz); error('Invalid bz shape.'); end
                % [nx,ny,nz] kpoint point density along each direction
                fid = fopen(fkpoints,'w');
                    fprintf(fid,'%s \n','KPOINTS');
                    fprintf(fid,'%s \n','0');
                    if all(bz==1); fprintf(fid,'%s \n','Gamma');
                    else;          fprintf(fid,'%s \n','Monkhorst-Pack'); end
                    fprintf(fid,' %i %i %i \n',ceil(bz));
                    fprintf(fid,' %i %i %i \n',[0 0 0]);
               fclose(fid);
            elseif isstring(bz)
                switch bz
                    case {'gamma','Gamma'}; am_dft.write_kpoints([1 1 1],fkpoints);
                    otherwise; error('Invalid bz flag.');
                end
            elseif isstruct(bz)
                if    isfield(bz,'natoms') % check if unit cell
                    % % still need to implement k-point density
                    % % calculates mesh for as isometric as possible k-point sampling
                    % % for VN: 10 = [ 5, 5, 5]; 20 = [ 9, 9, 9]; 30 = [13,13,13]; 40 = [17,17,17];
                    % ceil to nearest odd integer
                    kpt_ = @(bas,s) 2*ceil(am_lib.normc_(inv(bas).')*s/2)+1; n = kpt_(bz.bas,kdensity);
                    am_dft.write_kpoints(n,fkpoints);
                elseif isfield(bz,'nks') % check if brillouin zone
                    if isfield(bz,'w')
                        weights = bz.w;
                    else
                        weights = zeros(1,bz.nks); weights(:) = 1;
                    end
                    % explicit kpoints based
                    fid = fopen(fkpoints,'w');
                        fprintf(fid,'%s\n','KPOINTS');
                        fprintf(fid,'%i\n',bz.nks);
                        fprintf(fid,'%s\n',bz.units);
                        for i = 1:bz.nks
                            fprintf(fid,' %20.18f %20.18f %20.18f    %i\n',bz.k(:,i),weights(i));
                        end
                   fclose(fid);
                else
                    error('Invalid structure.');
                end
            end
        end

        function           write_qsub(hardware,account,nnodes,hrs,job,email)
            % defaults
            if nargin < 2; account = []; end
            if nargin < 3; nnodes = 1; end
            if nargin < 4; hrs = 12; end
            if nargin < 5; job = 'job'; end
            if nargin < 6; email = 'amei2@illinois.edu'; end

            fid = fopen('QSUB','w');
                switch hardware
                    case 'taub'
                        fprintf(fid,'#PBS -l nodes=%i:ppn=12\n',nnodes);
                        fprintf(fid,'#PBS -N %s \n',job);
                        fprintf(fid,'#PBS -l walltime=%i:00:00\n',hrs);
                        fprintf(fid,'#PBS -q cse \n');
                        fprintf(fid,'#PBS -V \n');
                        fprintf(fid,'module load openmpi/1.4-gcc\n');
                        fprintf(fid,'module load intel/11.1\n');
                        fprintf(fid,'cd \\${PBS_O_WORKDIR}\n');
                        fprintf(fid,'mpiexec /home/amei2/vasp.5.3.3 | tee OUTPUT \n');
                    case 'trio'
                        if isempty(account); error('Account is required for running on triolith'); end
                        fprintf(fid,'#!/bin/bash \n');
                        fprintf(fid,'#SBATCH -N %i \n',nnodes);
                        fprintf(fid,'#SBATCH -t %i:00:00 \n',hrs);
                        fprintf(fid,'#SBATCH -J %s \n',job);
                        fprintf(fid,'#SBATCH --exclusive \n');
                        fprintf(fid,'#SBATCH -A %s \n',account);
                        fprintf(fid,'#SBATCH --mail-user=%s \n',email);
                        fprintf(fid,'module load vasp/5.4.4-18Apr17 \n');
                        fprintf(fid,'module load python/2.7.12 \n');
                        fprintf(fid,'\n');
                        fprintf(fid,'# vasp_raman.py variables \n');
                        fprintf(fid,'export VASP_RAMAN_RUN=''mpprun vasp_std | tee OUTPUT'' \n');
                        fprintf(fid,'export VASP_RAMAN_PARAMS=''01_12_2_0.01'' \n');
                        fprintf(fid,'\n');
                        fprintf(fid,'#python ../_scripts/vasp_raman.py | tee RAMAN \n');
                        fprintf(fid,'mpprun vasp_std | tee OUTPUT \n');
                        fprintf(fid,'#mpprun vasp-gamma | tee OUTPUT \n');
                end
            fclose(fid);
        end

        function [x]  =    get_vasp(flag)
            import am_dft.*
            
            % check file is present
            flist = {'PROCAR','OUTCAR','vasprun.xml'};
            for i = 1:numel(flist)
                if     contains(flag,lower(flist{i}(1:6)))
                    if exist(flist{i},'file')~=2
                        fprintf('ERROR: Trying to read from %s when it does not exist',flist{i});
                        x = [];
                        return;
                    end
                end
            end
            
            switch flag
                % INCAR
                case 'incar:potim'
                    [~,x] = system('awk ''/POTIM/ { print $3 }'' INCAR');  x = sscanf(x,'%f');
                case 'incar:nedos'
                    [~,x] = system('awk ''/nedos/ { print $3 }'' INCAR');  x = sscanf(x,'%f');
                    
                % OSZICAR
                case 'oszicar:toten'
                    [~,x] = system('tail -n 1 OSZICAR | sed ''s/=/ = /g'' | awk ''{print $7}''');  x = sscanf(x,'%f');
                    
                % OUTCAR
            	case 'outcar:site_occupancy' % for DFT+U
                    [~,m] = system('awk ''/onsite density matrix/'' OUTCAR'); m = numel(strfind(m,'matrix'));
                    [~,x] = system('awk ''/ o = /{ print $3 }'' OUTCAR'); x = sscanf(x,'%f'); x = reshape(x,[],m);
                case 'outcar:toten'
                    [~,x] = system('awk ''/TOTEN/ { print $5 }'' OUTCAR');  x = sscanf(x,'%f');
                case 'outcar:entropy'
                    [~,x] = system('awk ''/EENTRO/ { print $5 }'' OUTCAR');  x = sscanf(x,'%f');
                case 'outcar:fermi'
                    [~,x] = system('awk ''/E-fermi/ { print $3 }'' OUTCAR');  x = sscanf(x,'%f');
                case 'outcar:pressure'
                    [~,x] = system('awk ''/pressure/ { print $4 }'' OUTCAR'); x = sscanf(x,'%f');
                case 'outcar:temperature'
                    [~,x] = system('awk ''/EKIN_LAT/ { print $6 }'' OUTCAR'); x = sscanf(x,'%f');
                case 'outcar:nedos'
                    [~,x] = system('awk ''/NEDOS/ { print $6 }'' OUTCAR'); x = sscanf(x,'%f');
                case 'outcar:nions'
                    [~,x] = system('awk ''/NIONS/ { print $12 }'' OUTCAR'); x = sscanf(x,'%f');
                case 'outcar:ispin'
                    [~,x] = system('awk ''/ISPIN/ { print $3 }'' OUTCAR'); x = sscanf(x,'%f');
                case 'outcar:nbands'
                    [~,x] = system('awk ''/NBANDS/ {print $15}'' OUTCAR'); x = sscanf(x,'%f');
                case 'outcar:kpoint'
                    [~,x] = system('awk ''/ plane waves: / { print $4 "  " $5 "  " $6}'' OUTCAR'); x = sscanf(x,'%f'); x = reshape(x,3,[]).';
                case 'outcar:gamma'
                    [~,x] = system('awk ''/GAMMA/ {print $4}'' OUTCAR'); x = sscanf(x,'%f');
                case 'outcar:niterations'
                    [~,x] = system('awk ''/Iteration/'' OUTCAR'); x = numel(strfind(x,'Iteration'));
                case 'outcar:magnetization:iteration'
                    [~,x] = system('awk ''/magnetization/'' OUTCAR | awk ''/electron/ { print $6 }'''); x = sscanf(x,'%f');
                case 'outcar:magnetization:ion'
                    natoms = get_vasp('vasprun:natoms');
                    [~,x] = system(sprintf('grep -A %i ''magnetization (x)'' OUTCAR | tail -n %i',natoms+3,natoms)); x = sscanf(x,'%f'); x = reshape(x,[],natoms).'; 
                    x = x(:,2:(end-1)); % remove index and total
                case 'outcar:phonon_energies' % obtained when ibrion=8 is specified
                    natoms = get_vasp('vasprun:natoms');
                    [~,x] = system(sprintf('grep 2PiTHz OUTCAR | head -n %i | sed ''s/=/ = /g'' | awk ''{print $10}'' ',3*natoms)); x = sscanf(x,'%f');
                case 'outcar:phonon_eigenvectors' % obtained with ibrion=8
                    natoms = get_vasp('vasprun:natoms');
                    [~,x] = system(sprintf('grep -A %i ''2PiTHz'' OUTCAR | grep -v X | grep -v 2PiTHz | awk ''{print $4 " " $5 " " $6}''',natoms+1)); 
                    x = sscanf(x,'%f'); x = x(1:(3*natoms).^2); x = reshape(x,3,natoms,3*natoms); % these ARE NOT divided by sqrt(m).
                case 'outcar:phonon_eigendisplacements' % obtained with ibrion=8 and iwrite = 3
                    natoms = get_vasp('vasprun:natoms');
                    [~,x] = system(sprintf('grep -A 100000 ''Eigenvectors after division by SQRT(mass)'' OUTCAR | grep -A %i ''2PiTHz'' | grep -v X | grep -v 2PiTHz | awk ''{print $4 " " $5 " " $6}''',natoms+1)); 
                    x = sscanf(x,'%f'); x = x(1:(3*natoms).^2); x = reshape(x,3,natoms,3*natoms); % these ARE divided by sqrt(m).
                case 'outcar:elastic' % obtained with ibrion=6 and iwrite = 3, the factor of 10 converts from kbar to GPa, the full tensor is given (not in voigt notation!)
                    [~,x] = system('grep -A 8 "TOTAL ELASTIC" OUTCAR | tail -n 6 | awk ''{print $2 " " $3 " " $4 " " $5 " " $6 " " $7}'' '); x = reshape(sscanf(x,'%f'),6,6).'/10; x = cij2cijkl(x);
                case 'outcar:total_charge' % obtained when iorbit=11 is specified
                    % x( natoms , {s,p,d,f} )
                    natoms = get_vasp('vasprun:natoms');
                    [~,x] = system(sprintf('grep -A %i ''total charge'' OUTCAR | tail -n %i | awk ''{print $2 " " $3 " " $4 " " $5}''',natoms+3,natoms)); x = sscanf(x,'%f'); x = reshape(x,[],natoms).';
                case 'outcar:real_df'
                    nedos = get_vasp('outcar:nedos');
                    [~,x] = system(sprintf('grep -A %i ''REAL DIELECTRIC FUNCTION'' OUTCAR | tail -n %i | awk ''{print $1 " " $2 " " $3 " " $4 " " $5 " " $6 " " $7}''',nedos+2,nedos)); x = sscanf(x,'%f'); x = reshape(x,[],nedos).';
                case 'outcar:imag_df'
                    nedos = get_vasp('outcar:nedos');
                    [~,x] = system(sprintf('grep -A %i ''IMAGINARY DIELECTRIC FUNCTION'' OUTCAR | tail -n %i | awk ''{print $1 " " $2 " " $3 " " $4 " " $5 " " $6 " " $7}''',nedos+2,nedos)); x = sscanf(x,'%f'); x = reshape(x,[],nedos).';
                case 'outcar:jdos'
                    nedos = get_vasp('outcar:nedos');
                    [~,x] = system(sprintf('grep -A %i ''JDOS'' OUTCAR | tail -n %i | awk ''{print $1 " " $2 " " $3 " " $4 " " $5 " " $6 " " $7}''',nedos+2,nedos)); x = sscanf(x,'%f'); x = reshape(x,[],nedos).';
                case 'outcar:plasma_frequency_inter'
                    [~,x] = system(sprintf('grep -A %i ''interband'' OUTCAR | tail -n %i | awk ''{print $1 " " $2 " " $3}''',3+2,3)); x = sscanf(x,'%f'); x = reshape(x,[],3).';
                % POTCAR
                case 'potcar:zval'  % valence on each atom in potcar
                    [~,x] = system('awk ''/ZVAL/{print $6}'' POTCAR'); x = sscanf(x,'%f');
                case 'potcar:vrhfin:configuration' % ionic configuration of each atom (orbitals + occupation)
                    [~,x] = system('grep ''VRHFIN'' POTCAR | sed -E ''s/([0-9])/\1,/g'' | sed -E ''s/ //g'' | sed -E ''s/,/ /g'' | sed -n -E ''s/^.*://p'' '); 
                    x = strtrim(strsplit(strtrim(x),'\n'));
                case 'potcar:vrhfin:orbital' % orbitals of each atom 
                    [~,x] = system('grep ''VRHFIN'' POTCAR | sed -E ''s/([0-9])//g'' | sed -E ''s/ //g'' | sed -n -E ''s/^.*://p'' '); 
                    x = strtrim(strsplit(strtrim(x),'\n'));
                    
                % PROCAR
                case 'procar:ispin'
                    [~,x] = system('awk ''/# of k-points:/'' PROCAR'); x = numel(strfind(x,'ions'));
                case 'procar:nkpts'
                    [~,x] = system('grep -m 1 "# of k-points:" PROCAR | awk ''{print $4}'''); x = sscanf(x,'%f');
                case 'procar:nbands'
                    [~,x] = system('grep -m 1 "# of k-points:" PROCAR | awk ''{print $8}'''); x = sscanf(x,'%f');
                case 'procar:natoms'
                    [~,x] = system('grep -m 1 "# of k-points:" PROCAR | awk ''{print $12}'''); x = sscanf(x,'%f');
                case 'procar:norbitals'
                    [~,x] = system('grep -m 1 py PROCAR'); x = numel(strsplit(strtrim(x))) - 2;
                case 'procar:orbital'
                    [~,x] = system('grep -m 1 py PROCAR'); x = strsplit(strtrim(x)); x=x(2:(end-1));
                case 'procar:kpoint'
                    [~,x] = system('grep weight PROCAR | sed ''s/-/ -/g'' | awk ''{print $5 " " $6 " " $7}'''); x = sscanf(x,'%f'); x = reshape(x,3,[]).';
                case 'procar:weight'
                    [~,x] = system('grep weight PROCAR | grep ''k-point'' | sed ''s/-/ -/g'' | awk ''{print $10}''' ); x = sscanf(x,'%f');
                case 'procar:energy'
                    nkpts     = get_vasp('procar:nkpts');
                    nbands    = get_vasp('procar:nbands');
                    nspins    = get_vasp('procar:ispin');
                    [~,x] = system('awk ''/energy/{print $5}'' PROCAR'); x = sscanf(x,'%f');
                    % E(nbands,nkpts,nspins)
                    x = reshape(x,nbands,nkpts,nspins);
                case 'procar:occupation'
                    nkpts     = get_vasp('procar:nkpts');
                    nbands    = get_vasp('procar:nbands');
                    nspins    = get_vasp('procar:ispin');
                    [~,x] = system('awk ''/energy/{print $8}'' PROCAR'); x = sscanf(x,'%f');
                    % E(nbands,nkpts,nspins)
                    x = reshape(x,nbands,nkpts,nspins);
                case 'procar:lmproj'
                    natoms    = get_vasp('procar:natoms');
                    norbitals = get_vasp('procar:norbitals');
                    nkpts     = get_vasp('procar:nkpts');
                    nbands    = get_vasp('procar:nbands');
                    nspins    = get_vasp('procar:ispin');
                    switch norbitals
                        case 1+3+5
                            [~,x]=system(sprintf('grep -A%i py PROCAR | grep -A%i tot | grep -v py | awk ''{print $2 " " $3 " " $4 " " $5 " " $6 " " $7 " " $8 " " $9 " " $10}''',natoms,natoms));
                        case 1+3+5+7
                            [~,x]=system(sprintf('grep -A%i py PROCAR | grep -A%i tot | grep -v py | awk ''{print $2 " " $3 " " $4 " " $5 " " $6 " " $7 " " $8 " " $9 " " $10 " " $11 " " $12 " " $13 " " $14 " " $15 " " $16 " " $17}''',natoms,natoms));
                    end
                    x = strrep(x,'--',''); x = sscanf(x,'%f'); x = reshape(x,norbitals,natoms,nbands,nkpts,nspins); % last is ispin
                    % lmproj(nspins,norbitals,nions,nbands,nkpts)
                    x = permute(x,[5,1,2,3,4]);
                case 'procar:phasefactor'
                    natoms    = get_vasp('procar:natoms');
                    norbitals = get_vasp('procar:norbitals');
                    nkpts     = get_vasp('procar:nkpts');
                    nbands    = get_vasp('procar:nbands');
                    nspins    = get_vasp('procar:ispin');
                    switch norbitals
                        case 1+3+5
                            [~,x]=system(sprintf('grep -v tot PROCAR | grep -A%i py | grep -v py | awk ''{print $2 " " $3 " " $4 " " $5 " " $6 " " $7 " " $8 " " $9 " " $10}''',2*natoms));
                        case 1+3+5+7
                            [~,x]=system(sprintf('grep -v tot PROCAR | grep -A%i py | grep -v py | awk ''{print $2 " " $3 " " $4 " " $5 " " $6 " " $7 " " $8 " " $9 " " $10 " " $11 " " $12 " " $13 " " $14 " " $15 " " $16 " " $17}''',2*natoms));
                    end
                    x = strrep(x,'--',''); x = strrep(x,'-',' -'); x = sscanf(x,'%f'); x = reshape(x,norbitals,2,natoms,nbands,nkpts,nspins); 
                    % add real and imaginary parts
                    x = permute( x(:,1,:,:,:,:) + 1i * x(:,2,:,:,:,:) ,[1,3,4,5,2] );
                    % lmproj(nspins,norbitals,nions,nbands,nkpts)
                    x = permute(x,[5,1,2,3,4]);
                case 'proout:nkpts'
                    [~,x] = system('grep -m 1 "# of k-points:" PROOUT.1 | awk ''{print $4}'''); x = sscanf(x,'%f');
                case 'proout:nbands'
                    [~,x] = system('grep -m 1 "# of k-points:" PROOUT.1 | awk ''{print $8}'''); x = sscanf(x,'%f');
                case 'proout:natoms'
                    [~,x] = system('grep -m 1 "# of k-points:" PROOUT.1 | awk ''{print $12}'''); x = sscanf(x,'%f');
                case 'proout:augmentation'
                    natoms    = get_vasp('proout:natoms');
                    norbitals = 9; % hard coded in sphpro.F
                    nkpts     = get_vasp('proout:nkpts');
                    nbands    = get_vasp('proout:nbands');
                    nspins    = get_vasp('procar:ispin');
                    [~,x]=system(sprintf('grep -A%i augmentation PROOUT.1 | grep -v augmentation | awk ''{print $1 " " $2 " " $3 " " $4 " " $5 " " $6 " " $7 " " $8 " " $9}''',natoms*nbands*nkpts));
                    %
                    x = strrep(x,'--',''); x = strrep(x,'-',' -'); x = sscanf(x,'%f'); x = reshape(x,norbitals,natoms,nbands,nkpts,nspins); 
                    % augmentation(nspins,norbitals,natoms,nbands,nkpts)
                    x = permute(x,[5,1,2,3,4]);
                case 'proout:phase'
                    natoms    = get_vasp('proout:natoms');
                    norbitals = 9; % hard coded in sphpro.F
                    nkpts     = get_vasp('proout:nkpts');
                    nbands    = get_vasp('proout:nbands');
                    nspins    = 1; % multiple spins not yet implemented
                    [~,x]=system(sprintf('grep -B%i augmentation PROOUT.1 | grep -v augmentation | awk ''{print $1 " " $2 " " $3 " " $4 " " $5 " " $6 " " $7 " " $8 " " $9}''',2*nkpts*nbands*natoms));
                    % 
                    x = strrep(x,'--',''); x = strrep(x,'-',' -'); x = sscanf(x,'%f'); x = reshape(x,2,norbitals,nbands,nkpts,natoms,nspins); 
                    % add real and imaginary parts
                    x = x(1,:,:,:,:,:) + 1i * x(2,:,:,:,:,:);
                    % augmentation(nspins,norbitals,nions,nbands,nkpts)
                    x = permute(x,[1,2,5,3,4]);
                % VASPRUN
                case 'vasprun:scf_energy'
                    [~,x] = system('grep ''e_fr_energy'' vasprun.xml | sed ''s/\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*/ 1e8 /g'' | awk ''{print $3}'''); x = sscanf(x,'%f');
                case 'vasprun:natoms'
                    [~,x] = system('grep ''<atoms>'' vasprun.xml | awk ''{print $2}'''); x = sscanf(x,'%f');
                otherwise 
                    error('grep command not found');
            end
            if isempty(x); x = NaN; end
            

            % convert voigt elastic constants to full tensor
            function[CC] = cij2cijkl(C)
                CC=zeros(3,3,3,3);
                CC(1,1,1,1)=C(1,1);CC(2,2,2,2)=C(2,2);CC(3,3,3,3)=C(3,3);CC(2,3,2,3)=C(4,4);CC(3,2,3,2)=CC(2,3,2,3);CC(2,3,3,2)=CC(2,3,2,3);
                CC(3,2,2,3)=CC(2,3,2,3);CC(1,3,1,3)=C(5,5);CC(3,1,1,3)=CC(1,3,1,3);CC(1,3,3,1)=CC(1,3,1,3);CC(3,1,3,1)=CC(1,3,1,3);
                CC(1,1,2,2)=C(1,2);CC(2,2,1,1)=CC(1,1,2,2);CC(1,1,3,3)=C(1,3);CC(3,3,1,1)=CC(1,1,3,3);CC(1,1,2,3)=C(1,4);CC(1,1,3,2)=CC(1,1,2,3);
                CC(2,3,1,1)=CC(1,1,2,3);CC(3,2,1,1)=CC(1,1,2,3);CC(1,1,1,3)=C(1,5);CC(1,1,3,1)=CC(1,1,1,3);CC(1,3,1,1)=CC(1,1,1,3);CC(3,1,1,1)=CC(1,1,1,3);
                CC(1,1,1,2)=C(1,6);CC(1,1,2,1)=CC(1,1,1,2);CC(1,2,1,1)=CC(1,1,1,2);CC(2,1,1,1)=CC(1,1,1,2);CC(2,2,3,3)=C(2,3);CC(3,3,2,2)=CC(2,2,3,3);
                CC(2,2,2,3)=C(2,4);CC(2,2,3,2)=CC(2,2,2,3);CC(2,3,2,2)=CC(2,2,2,3);CC(3,2,2,2)=CC(2,2,2,3);CC(2,2,1,3)=C(2,5);CC(2,2,3,1)=CC(2,2,1,3);
                CC(1,3,2,2)=CC(2,2,1,3);CC(3,1,2,2)=CC(2,2,1,3);CC(2,2,1,2)=C(2,6);CC(2,2,2,1)=CC(2,2,1,2);CC(1,2,2,2)=CC(2,2,1,2);CC(2,1,2,2)=CC(2,2,1,2);
                CC(3,3,2,3)=C(3,4);CC(3,3,3,2)=CC(3,3,2,3);CC(2,3,3,3)=CC(3,3,2,3);CC(3,2,3,3)=CC(3,3,2,3);CC(3,3,1,3)=C(3,5);CC(3,3,3,1)=CC(3,3,1,3);
                CC(1,3,3,3)=CC(3,3,1,3);CC(3,1,3,3)=CC(3,3,1,3);CC(3,3,1,2)=C(3,6);CC(3,3,2,1)=CC(3,3,1,2);CC(1,2,3,3)=CC(3,3,1,2);CC(2,1,3,3)=CC(3,3,1,2);
                CC(2,3,1,3)=C(4,5);CC(3,2,1,3)=CC(2,3,1,3);CC(1,3,3,2)=CC(2,3,1,3);CC(1,3,2,3)=CC(2,3,1,3);CC(2,3,3,1)=CC(2,3,1,3);CC(3,2,3,1)=CC(2,3,1,3);
                CC(3,1,2,3)=CC(2,3,1,3);CC(3,1,3,2)=CC(2,3,1,3);CC(2,3,1,2)=C(4,6);CC(3,2,1,2)=CC(2,3,1,2);CC(1,2,2,3)=CC(2,3,1,2);CC(1,2,3,2)=CC(2,3,1,2);
                CC(2,3,2,1)=CC(2,3,1,2);CC(3,2,2,1)=CC(2,3,1,2);CC(2,1,2,3)=CC(2,3,1,2);CC(2,1,3,2)=CC(2,3,1,2);CC(1,3,1,2)=C(5,6);CC(3,1,1,2)=CC(1,3,1,2);
                CC(1,2,1,3)=CC(1,3,1,2);CC(1,2,3,1)=CC(1,3,1,2);CC(1,3,2,1)=CC(1,3,1,2);CC(3,1,2,1)=CC(1,3,1,2);CC(2,1,1,3)=CC(1,3,1,2);
                CC(2,1,3,1)=CC(1,3,1,2);CC(1,2,1,2)=C(6,6);CC(2,1,1,2)=CC(1,2,1,2);CC(1,2,2,1)=CC(1,2,1,2);CC(2,1,2,1)=CC(1,2,1,2);
            end
        end
            
        function           vasp_fix_input(flags)
            
            import am_dft.* am_lib.*
            
            % no flag
            if nargin < 1; flags=''; end
            
            % defaults
            if contains(flags,'dryrun'); isdryrun = true; else; isdryrun=false; end

            % PRIMITIVE BASIS UNREASONABLY LARGE
            uc = load_poscar('POSCAR'); d = uc.natoms/abs(det(uc.bas*0.1));
            if d<0.001 
                fprintf('DIR: %s \n',pwd); 
                fprintf('       ERROR: at least of the lattice vectors is too long. \n');
                fprintf('       bas: %10.5f %10.5f %10.5f \n',uc.bas(:,1));
                fprintf('            %10.5f %10.5f %10.5f \n',uc.bas(:,2));
                fprintf('            %10.5f %10.5f %10.5f \n',uc.bas(:,3));
                target_density = 50; % [atoms/nm^3 ~ 5E22 atoms/cm3]
                uc.bas = (d/target_density)^(1/3)*uc.bas;
                fprintf('       FIX: scaling lattice vectors to get a reasonable density of ~ %i atoms/nm^3 \n',target_density)
                fprintf('       bas: %10.5f %10.5f %10.5f \n',uc.bas(:,1));
                fprintf('            %10.5f %10.5f %10.5f \n',uc.bas(:,2));
                fprintf('            %10.5f %10.5f %10.5f \n',uc.bas(:,3));
                if ~isdryrun
                    write_poscar(uc);
                end
            end

            % SCF DIVERGED
            b = get_vasp('vasprun:scf_energy');
            if ~isempty(b) && any(abs(b)>1E7)
                fprintf('DIR: %s \n',pwd);
                fprintf('       ERROR: electronic cycle divreged (d eps = %g). \n',b(end));
                % fix is to try a lower AMIX and set SCF algo to Normal
                incar = load_incar('INCAR');
                % reduce amix by 1/2
                if isfield(incar,'amix') && ~isempty(incar.amix); amix = str2double(incar.amix); else; amix = 0.4; end
                amix = amix/2;
                incar = generate_incar(incar,{'algo','Normal','amix',num2str(amix)});
                % print the fix
                fprintf('       FIX: set algo = normal and reduce amix to %0.2f\n',amix)
                if ~isdryrun
                    write_incar(incar);
                end
                return;
            end
            
            % NOT ENOUGH BANDS
            if grepcontains_('Your highest band is occupied at some k-points!','output')
                fprintf('DIR: %s \n',pwd); 
                nbands = get_vasp('outcar:nbands');
                fprintf('       ERROR: not enough bands (nbands = %g). \n',nbands);
                % read incar
                incar = load_incar('INCAR');
                if isfield(incar,'ncore'); ncore = str2double(incar.ncore); else; ncore = 1; end
                % update nbands
                nbands = ceil(nbands*1.2); nbands = ncore-mod(nbands,ncore)+nbands;
                fprintf('       FIX: increase nbands by ~1.2x to %i\n',nbands)
                % fix is to try a lower AMIX and set SCF algo to Normal
                incar = generate_incar(incar,{'nbands',nbands});
                if ~isdryrun
                    write_incar(incar);
                end
            end

            % NONHERMITIAN MATRIX
            if grepcontains_('WARNING: Sub-Space-Matrix is not hermitian!','output')
                fprintf('DIR: %s \n',pwd);
                fprintf('       ERROR: Sub-Space-Matrix is not hermitian. \n');
                fprintf('       This is usually caused by some other problem. Investigate the output manually. \n')
            end

            
            function b = grepcontains_(str,file)
                [~,b]=system(sprintf('grep ''%s'' %s',str,file));
                b = ~isempty(b);
            end
        end

        % vasp automation
        
        function           dfpt2frozen(pc,incar,d,hw,modelist,dmax,ndisps,flag)
            % clear;clc;odir=pwd;
            % % select a mode, maximum displacement, number of displacements
            % modelist = [1,2,4,10,16,22]; dmax = 0.01; ndisps = 15; flag = 'A';
            % % load phonon eigenvectors and unit cell
            % bdir='./DFPT';cd(bdir);
            %     [~,pc] = am_cell.load_cell('poscar','POSCAR');
            %     incar  = am_dft.load_incar('INCAR');
            %     d      = am_dft.get_vasp('outcar:phonon_eigendisplacements'); % obtained with iwrite = 3
            %     % V      = am_dft.get_vasp('outcar:phonon_eigenvectors');
            %     % d      = V./sqrt(pc.mass(pc.species));
            %     hw     = am_dft.get_vasp('outcar:phonon_energies');
            % cd(odir);
            % dfpt2frozen(pc,incar,d,hw,modelist,dmax,ndisps,flag);

            % set origin directory
            odir=pwd;
            % check modelist
            modelist = modelist(modelist<=size(d,3));
            % get scale factors setting maximum displacement
            scale_ = @(d,imode) max(am_lib.normc_(d(:,:,imode)));
            % create directories
            wdir_ = @(imode,idisp) sprintf('./%.3i/%.3i',imode,idisp);
            % setup or analysis?
            switch flag
                case {'S','setup'}
                    % modify incar
                    incar.ibrion = '-1'; incar.npar  =  ''; 
                    incar.kpar   = '8';  incar.ncore = '8';
                    % set the maximum displacement to between [-0.01,0.01]; % [nm]
                    alpha = linspace(-dmax,dmax,ndisps).'/scale_(d,imode);
                    % loop
                    for i = 1:numel(modelist); imode = modelist(i);
                    for idisp = 1:numel(alpha)
                        wdir=wdir_(imode,idisp); mkdir(wdir); cd(wdir);
                            % write files
                            % POSCAR
                                sc = pc; sc.tau = sc.tau + sc.bas\(alpha(idisp)*d(:,:,imode));
                                am_dft.write_poscar(sc,'POSCAR')
                            % INCAR
                                am_dft.write_incar(incar);
                            % KPOINTS, POTCAR
                            for f = {'KPOINTS','POTCAR'}
                                copyfile([bdir,'/',f{:}],f{:});
                            end
                        cd(odir);
                    end
                    end
                case {'A','analysis'}
                    figure(1);clc;clf;set(gcf,'color','w');hold on;
                    for i = 1:numel(modelist); imode = modelist(i);
                        % grep oszicar for final energy
                        E = NaN(ndisps,1);
                        for idisp = 1:ndisps
                            wdir=wdir_(imode,idisp); if exist(wdir,'dir')~=7; continue; end; cd(wdir);
                                E(idisp) = am_dft.get_vasp('oszicar:toten');
                            cd(odir);
                        end
                        % set the maximum displacement to between [-0.01,0.01]; % [nm]
                        alpha = linspace(-dmax,dmax,ndisps).'./scale_(d,imode);
                        % plot and fit results
                        if sum(~isnan(E))>2
                            scale = max(am_lib.normc_(d(:,:,imode))); ex_ = ~isnan(E);
                            % U [eV] = 1/2 * u^2 [nm/sqrt(amu)] * amu;
                            plot(alpha(ex_),E(ex_),'o-');
                            ft_ = fit(alpha(ex_),E(ex_),'poly2');
                            plot(alpha,E,'o-',alpha,feval(ft_,alpha),'--');

                            phonon_energy = sqrt( 2 * ft_.p1 * am_dft.units_meV * 2*pi );
                            [ phonon_energy, hw(imode), phonon_energy./hw(imode) ]
                        end
                    end
            end
        end
        
    end
    
    
    % core

    methods (Static)
        
        % test suite
        
        function test
            
            import am_dft.* am_lib.*
            
            % test space group memory is up to date
            for i = 1:237
                [S1] = generate_sg(i,true);
                [S2] = generate_sg(i,false);
                criteria_S(i) = all_(eq_(S1,S2));
            end
            test_(all(criteria_S),'sg memory check',sprintf('failed to match symmetry with memory for sg:%s',sprintf(' %i',find(~criteria_S))));
            
            % test point group memory is up to date
            for i = 1:32
                [R1,W1] = generate_pg(i,true);
                [R2,W2] = generate_pg(i,false);
                criteria_R(i) = all_(eq_(R1,R2));
                criteria_W(i) = all_(eq_(W1,W2));
            end
            test_(all(criteria_R),'pg memory check',sprintf('failed to match symmetry with memory for pg:%s',sprintf(' %i',find(~criteria_R))));
            test_(all(criteria_W),'dg memory check',sprintf('failed to match symmetry with memory for dg:%s',sprintf(' %i',find(~criteria_W))));
            
            % generate and test point groups
            for i = 1:32
                % recover pg_code
                R = generate_pg(i,true); pg_id(i) = identify_pointgroup(R); 
                % make sure number of classes match number of irreps
                CT = get_character_table( get_irreducible_representations(R) ); [nirreps(i),nclasses(i)] = size(CT);
            end
            criteria = pg_id==[1:32];
            test_(all(criteria),'pg generation',sprintf('failed to generate pointgroups:%s',sprintf(' %i',find(~criteria))));
            criteria = nirreps==nclasses;
            test_(all(criteria),'pg classes+irreps',sprintf('failed to produce an equal number of classes and irreps:%s',sprintf(' %i',find(~criteria))));
            
            % test winger j=1 rotation matrices
            R = cat(3,generate_pg(32,false),generate_pg(27,false));
            W = get_wigner(1,R,'tesseral');
            criteria = squeeze(all(all(eq_(W,R),1),2));
            test_(all(criteria),'wigner SO(3)',sprintf('failed to generate consistent SO(3) rotations:%s',sprintf(' %i',find(~criteria))));
            
            % test double groups
            for i = 1:32
                % get single-valued representation
                [R,W] = generate_pg(i,false);
                % compare number of symmetries
                nRs(i) = size(R,3);
                nWs(i) = size(complete_group(W),3);
                % compare 
                CT = get_character_table( get_irreducible_representations(W) ); [nirreps(i),nclasses(i)] = size(CT);
            end
            criteria = 2*nRs==nWs;
            test_(all(criteria),'dg generation',sprintf('failed to generate correct number of symmetries for dg:%s',sprintf(' %i',find(~criteria))));
            criteria = nirreps==nclasses;
            test_(all(criteria),'dg classes+irreps',sprintf('failed to produce an equal number of classes and irreps for dg:%s',sprintf(' %i',find(~criteria))));
            
            function test_(logical,test_name,fail_msg)
                if logical
                    fprintf('      %s: pass\n',test_name);
                else
                    fprintf('      %s: %s\n',test_name,fail_msg);
                end 
            end
        end
        


        % symmetry

%         function [sg,pg,dg]   = get_groups(pc, tol)
%
%             % getting point groups
%             [T,H,S,R]    = get_symmetries(pc, tol)
%
%             fprintf(' ... getting group properties')
%
%             % get multiplication table
%             [MT,E,I] = get_multiplication_table(S);
%             if all(all(MT-MT.'==0))
%                 fprintf(', abelian' );
%             else
%                 fprintf(', non-abelian' );
%             end
%
%             % get generators
%             G = get_generators(MT);
%             fprintf(', %i generators (', numel(G)); fprintf(' %i', G); fprintf(')');
%
%             % character table
%
%         end


        function [MT,E,I]     = get_multiplication_table(S,tol)
            % [MT,E,I] = get_multiplication_table(S)
            % get multiplication table: S(:,:,i)*S(:,:,j) = S(:,:,MT(i,j))
            
            import am_lib.* am_dft.*

            if nargin<2; tol=am_dft.tiny; end
            
            s = size(S); if numel(s)<3; s(3) = 1; end
            
            if iscell(S)  % seitz operator combined with permutation (represented as a two-part cell)
                s = size(S{1}); if numel(s)<3; s(3) = 1; end
                if     s(1)==4 && s(2)==4 && all_(eq_(S{1}(4,1:4,:), [0,0,0,1], tol)) % seitz operator (applies mod to translational components)
                    md_ = @(X) [X(1:12,:);mod_(X(13:15,:),tol);X(16:end,:)];
                    rs_ = @(X) md_(reshape(X,s(1)*s(2),[])); nsyms = s(3);

                    ref = [rs_(S{1});reshape(S{2},size(S{2},1),[])];
                    opr = [rs_(matmulp_(S{1},permute(S{1},[1,2,4,3]))); reshape(operm_(S{2},S{2}),size(S{2},1),[])];

                    MT = reshape( member_( opr , ref, tol), s(3), s(3)); 
                elseif s(1)==3 && s(2)==3                                          % point operator
                    rs_ = @(X) reshape(X,s(1)*s(2),[]); nsyms = s(3);

                    ref = [rs_(S{1});reshape(S{2},size(S{2},1),[])];
                    opr = [rs_(matmulp_(S{1},permute(S{1},[1,2,4,3]))); reshape(operm_(S{2},S{2}),size(S{2},1),[])];

                    MT = reshape( member_( opr , ref, tol), s(3), s(3));
                end
            elseif s(1)==4 && s(2)==4 && all_(eq_(S(4,1:4,:), [0,0,0,1], tol))     % seitz operator (applies mod to translational components)
                md_ = @(X) [X(1:12,:);mod_(X(13:15,:),tol);X(16:end,:)];
                rs_ = @(X) md_(reshape(X,s(1)*s(2),[])); nsyms = s(3);
                
                MT = reshape( member_( rs_(matmulp_(S,permute(S,[1,2,4,3]))), rs_(S), tol) , nsyms, nsyms);
            else                                                                   % point operator
                rs_ = @(X) reshape(X,s(1)*s(2),[]); nsyms = s(3); 

                MT = reshape( member_( rs_(matmulp_(S,permute(S,[1,2,4,3]))), rs_(S), tol ) , nsyms, nsyms);
            end

            if any(sum(MT,1)~=sum([1:nsyms])) || any(sum(MT,2)~=sum([1:nsyms]))
                error('MT is incorrect. Check for mistakes in the symmetry and ensure that the symmetry is in the primitive basis');
            end

            % find identity
            if nargout>1; E = find(all(MT==[1:nsyms].',1)); end

            % get inverse indicies
            if nargout>2; I = [MT==E]*[1:nsyms].'; end
        end
        
        function [MT]         = relabel_multiplication_table(MT,fwd)
            % relabel multiplication table
            rev(fwd) = [1:size(MT,1)]; MT = rev(MT(fwd,fwd));
        end
        
        function [subgroup,grph] = get_subgroups(MT)
            % clear;clc; import am_dft.* am_lib.*
            % [R,~]=generate_pg('o',false);
            % % get multiplication table
            % MT=am_dft.get_multiplication_table(R);
            % % get subgroups
            % [subgroup,structure] = get_subgroups(MT);
            % % identify the subgroup
            % for i = 1:size(subgroup,2); pg_id(i) = identify_pointgroup(R(:,:,logical(subgroup(:,i)))); end
            % % plot the group/subgroup hierarchy
            % h=plot(digraph(structure,'OmitSelfLoops'),'Layout','layered','NodeLabel',decode_pg(pg_id));

            import am_lib.* am_dft.*
            
            % set maximum generators (3 should be enough for point groups)
            maxgens = 3;

            % gen number of symmerties
            nsyms = size(MT,1);

            % preallocate space
            maxmem = 0; for ngens = 1:maxgens; maxmem = maxmem+nchoosek(nsyms,ngens); end
            subgroup = false(nsyms,maxmem);

            % [gen_id] = get_generators(MT)
            flatten_ = @(x) x(:); k=0;
            for ngens = 1:maxgens
                % generate a list of all possible generators (ignores generator order for speed)
                genlist = nchoosek_(nsyms,ngens);
                % loop over the list one set of generators at a time, checking whether the entire multiplication table can be generated
                for i = 1:size(genlist,2)
                    % initialize, include identity
                    x = genlist(:,i);
                    % expand multiplication table
                    while true
                        u = unique([flatten_(MT(x,x));x],'stable');
                        switch numel(u)
                            case {nsyms, numel(x)}
                                % exausted all possibiltiies: a subgroup or the entire group has been generated
                                k=k+1; 
                                subgroup(u,k) = true;
                                break;
                            otherwise
                                % update
                                x = u;
                        end
                    end
                end
            end

            % get unique subgroups
            subgroup = uniquec_(single(subgroup)); % nsubs = size(subgroup,2);

            % sort subgroups by order
            subgroup = subgroup(:, rankc_(sum(subgroup,1)) );

            % structure(:,k) indicates subgroups contained within group k
            grph = eq_((subgroup./normc_(subgroup).^2).'*subgroup,1);

        end
        
        function [G,u]        = get_generators(MT)
            % [G,u] = get_generators(MT)
            % u reindexes symmetries based on the generators
            import am_lib.*
            % [gen_id] = get_generators(MT)
            nsyms = size(MT,1); flatten_ = @(x) x(:);
            for ngens = 1:nsyms
                % generate a list of all possible generators (ignores generator order for speed)
                genlist = nchoosek_(nsyms,ngens);
                % loop over the list one set of generators at a time, checking whether the entire multiplication table can be generated
                for i = 1:size(genlist,2)
                    % initialize, include identity
                    x = genlist(:,i);
                    % expand multiplication table
                    while true
                        u = unique([flatten_(MT(x,x));x],'stable');
                        switch numel(u)
                            case nsyms
                                % found!
                                G = genlist(:,i);
                                return;
                                %
                            case numel(x)
                                % exausted all possibiltiies
                                break;
                            otherwise
                                % update
                                x = u;
                        end
                    end
                end
            end
        end

        function [Glist,ulist]= get_generators_list(MT,ngens)
            % returns all possible sets of generators with ngens elements
            import am_lib.*
            % [gen_id] = get_generators(MT)
            nsyms = size(MT,1); flatten_ = @(x) x(:);
            % generate a list of all possible generators (incldues order of generator for completness)
            Glist = perm_norep_(nsyms,ngens); ngenlists = size(Glist,2);
            ex_ = false(1,ngenlists); ulist = zeros(nsyms,ngenlists);
            % loop over the list one set of generators at a time, checking whether the entire multiplication table can be generated
            for i = 1:ngenlists
                % initialize
                x = Glist(:,i);
                % expand multiplication table
                while true
                    u = unique([flatten_(MT(x,x));x],'stable');
                    switch numel(u)
                        case nsyms
                            % found!
                            ex_(i) = true;
                            ulist(:,i) = u;
                            break;
                        case numel(x)
                            ex_(i) = false;
                            % exausted all possibiltiies
                            break;
                        otherwise
                            % update
                            x = u;
                    end
                end
            end
            Glist = Glist(:,ex_); ulist = ulist(:,ex_);
        end

        function                plot_cayley_graph(MT)
            % clear;clc; import am_dft.* am_lib.*
            % [R,~]=generate_pg('o_h');
            % % get multiplication table
            % MT=am_dft.get_multiplication_table(R);
            % % plot cayley graph
            % plot_cayley_graph(MT)
            %
            import am_dft.* am_lib.*
            [G] = get_generators(MT);
            i=[]; j=[]; x = MT(1,:);
            for k = 1:numel(G)
                i = [i;flatten_(x)];
                j = [j;flatten_(MT(G(k),x))];
            end
            plot(digraph(i,j,'OmitSelfLoops'),'layout','subspace3');
        end

        function [RR]         = get_regular_representation(S)
            % get regular rep G by putting identity along diagonal of multiplciation table
            [MT,~,I] = get_multiplication_table(S); nSs = size(MT,2);
            RR = double(accessc_(MT,I)==permute([1:nSs],[1,3,2]));
        end
        
        function [IR]         = get_irreducible_representations(S)
            % s2c = identifies the class to which symmetries belong
            % CT = character table
            % irreps
            %
            % MATHEMATICS OF COMPUTATION, VOLUME 24, NUMBER 111, JULY, 1970
            % Computing Irreducible Representations of Groups
            % By John D. Dixon

            import am_lib.* am_dft.*
            
            % get regular rep G by putting identity along diagonal of multiplciation table
            [MT,~,I] = get_multiplication_table(S); nSs = size(MT,2);
            RR = double(accessc_(MT,I)==permute([1:nSs],[1,3,2]));
            
            % initialize decomposition loop
            [~,~,nSs] = size(RR); U = eye(nSs); inds = ones(nSs,1); ninds = 1;

            % loop until irreps are fully decomposed
            while true
                % loop over cycle structures
                for j = 1:max(inds)
                    ex_ = inds==j; 
                    H = dixon_decomposition_( RR(ex_,ex_,:) ); [Vp,~] = eig(H,'vector');
                    for ig = 1:nSs
                        RR(ex_,ex_,ig) = Vp \ RR(ex_,ex_,ig) * Vp;
                    end
                    U(:,ex_) = (U(:,ex_)*Vp);
                end
                A = merge_(sum(abs(RR),3)>am_lib.eps); inds = [1:size(A,1)]*A;
                if ninds == max(inds); break; else; ninds = max(inds); end
            end

            % get irreducible representations
            nirreps = max(inds); IR = cell(1,nirreps);
            for i = 1:nirreps
                IR{i} = RR(inds==i,inds==i,1:nSs);
            end
            
            % sort the irreps by size
            IR = IR(rankc_(cellfun(@(x)size(x,1),IR)));
            
            function H = dixon_decomposition_(rr)
                nbases=size(rr,1);
                nsyms =size(rr,3);

                for r = 1:nbases
                for s = 1:nbases
                    Hrs = zeros(nbases);
                    if     r==s
                        Hrs(r,s) = 1;
                    elseif r>s
                        Hrs(r,s) = 1;
                        Hrs(s,r) = 1;
                    elseif r<s
                        Hrs(r,s) = sqrt(-1);
                        Hrs(s,r) =-sqrt(-1);
                    end

                    H(1:nbases,1:nbases) = 0;
                    for q = 1:nsyms
                        H = H + rr(:,:,q) \ Hrs * rr(:,:,q);
                    end
                    H = H / nsyms;

                    if any(abs( H(1,1)*eye(nbases)-H ) >am_lib.eps); return; end
                end
                end
                H = eye(nbases);
            end
        end

        function [CT,cc_id,ir_id] = get_character_table(IR)
            % R=generate_pg(32,true);
            % IR = get_irreducible_representations(R);
            % [CT,cc_id,ir_id] = get_character_table(IR);
            % print_character_table(R,CT,cc_id)
            
            import am_lib.*
            
            %
            nirreps = numel(IR);
            nsyms   = size(IR{1},3);
            
            % get character table
            CT = zeros(nirreps,nsyms);
            for i = 1:nsyms; for j = 1:nirreps
                CT(j,i) = trace(IR{j}(:,:,i));
            end; end
        
            % correct numerical error
            CT = wdv_(CT);

            % get irreducible irreps
            [CT,~,cc_id]=uniquec_(CT);

            % get irreducible classes
            [CT,~,ir_id]=uniquec_(CT.'); CT=CT.';
            
            % check 
            if ~isdiag_(CT'*CT)
                sym(CT)
                error('Character table is not orthogonal')
            end
        end

        function                print_character_table(R,CT,cc_id)
            % clear;clc
            % R=generate_pg(32,true);
            % IR = get_irreducible_representations(R);
            % [CT,cc_id,ir_id] = get_character_table(IR);
            % print_character_table(R,CT,cc_id)

            
            %     clear;clc; import am_dft.* am_lib.*
            %     [R,W]=generate_pg('o_h',false);
            %     % d = permute(det_(R),[3,1,2]);
            %     % j=1/2; [W] = get_wigner(j,d,'spherical');
            % 
            %     % (1i*eq_(d,-1) + eq_(d,1));
            %     % W = cat(3,1i*W,-1i*W);
            %     % W = cat(3,kron_(R,eye(2)),kron_(R,[0 1; 1 0]));
            %     % W = complete_group(W);
            %     % MT = am_dft.get_multiplication_table(W);
            %     % MTk = am_dft.get_multiplication_table(K);
            %     IR = get_irreducible_representations(W);
            %     [CT,cc_id,ir_id] = get_character_table( IR );
            %     % sym(CT)
            %     print_character_table(W,CT,cc_id)
            
            import am_dft.* am_lib.*
            
            % identify the prototypical symmetries
            [~,cc_rep]=unique(cc_id,'stable');

            % name of prototypical symmetries
            ps_name_long = get_long_ps_name(R);
                
            % get number of classes and number of elements in each class
            nclasses = numel(cc_rep); nelems = sum(cc_id(cc_rep)==cc_id',2).'; nsyms = sum(nelems);
            
            % calculate orbital characters
            [chi,J] = get_orbital_characters(R(:,:,cc_rep));

            % print character table% class label 
                print_table_(nclasses,wdv_(CT),chi,ps_name_long,nelems,cc_id,cc_rep)
            
            % DECOMPOSITIONS
                decomp_ = @(chi) wdv_( chi * (CT.*nelems).' / nsyms );
                
            % IRREP x IRREP decomposition
                decomp = decomp_( CT.^2 ); decomp = wdv_(decomp); nirreps = size(CT,1);
                for i = 1:nirreps; row_name{i} = sprintf('irr # %g',i); end
                for i = 1:nirreps; col_name{i} = sprintf('irr # %g',i); end
                print_table_decomposition_('irrep x irrep', row_name, col_name, decomp);
                
            % ORBITAL decomposition
                decomp = decomp_( chi ); decomp = wdv_(decomp);
                for i = 1:nirreps; row_name{i} = sprintf('irr # %g',i); end
                clear row_name; for i = 1:numel(J); row_name{i} = sprintf('%s',sym(J(i))); end
                clear col_name; for i = 1:nirreps;  col_name{i} = sprintf('irr # %g',i); end
                print_table_decomposition_('orbitals', row_name, col_name, decomp);
                
                
            function [chi,J] = get_orbital_characters(W)
                % identify double groups symmetries (negative ps_id = double valued)
                [ps_id,~,d] = am_dft.identify_point_symmetries(W);
                % get rotation angle
                alpha = am_lib.R_angle_(permute(d,[3,1,2]).*W);
                % compute character
                character_ = @(l,alpha) sin((l+1/2).*(alpha+1E-12))./sin((alpha+1E-12)./2); 
                chi_l = [0.0:5.0]; chi_O3 = character_(chi_l(:), alpha) .* (am_lib.eq_(d,-1).*1i + 1*am_lib.eq_(d,1));% .* (d).^(chi_l(:));% .* sign(ps_id);
                chi_j = [0.5:5.5]; chi_U2 = character_(chi_j(:), alpha) .* (d).^(chi_j(:)) .* sign(ps_id);
                J = [chi_l(:);chi_j(:)]; chi = [chi_O3;chi_U2]; chi = am_lib.wdv_(chi);
            end
            
            function print_table_(nclasses,CT,chi,ps_name_long,nelems,cc_id,cc_rep)
                fprintf('      '); for j = 1:nclasses; fprintf('%12s',['#',num2str(j)]); end; fprintf('\n');
                print_bar(nclasses);
                % print class name 
                fprintf('      '); for j = 1:nclasses; fprintf('%12s',ps_name_long{cc_rep(j)}); end; fprintf('\n');
                % print class elements
                fprintf('      '); for j = 1:nclasses; fprintf('%12i',nelems(j)); end; fprintf('\n');
                print_bar(nclasses);
                % print character table
                for l = 1:size(CT,1); fprintf('      '); fprintf('%12s',sym(CT(l,:))); fprintf('\n'); end
                % print SO(3) and SU(2) characters single and double groups
                print_bar(nclasses);
                for l = 1:size(chi,1); fprintf('      '); fprintf('%12s',sym(chi(l,:))); fprintf('\n'); end
                print_bar(nclasses);

                % symmetries in each class
                    fprintf('      '); fprintf('Symmetries in each classes:\n')
                    for j = 1:nclasses
                        fprintf('%12s:%s\n',['#',num2str(j)], sprintf(' %-12s',ps_name_long{cc_id==j}));
                    end
            end
            
            function print_table_decomposition_(table_name,row_name,col_name,decomp)
                [nrows,~] = size(decomp); 
                fprintf('\n');
                fprintf('      '); fprintf('%12s\n',table_name); 
                fprintf('      '); fprintf('%12s\n','==========='); 

                fprintf('      '); fprintf('%12s',''); fprintf('%12s',col_name{:}); fprintf('\n');
                print_bar(numel(col_name)+1)
                % fprintf('      '); fprintf('%12s',repmat(' -----------',1,ncols+1)); fprintf('\n');
                for k = 1:nrows
                    fprintf('        '); fprintf('%-10s',row_name{k}); fprintf('%12g',decomp(k,:)); fprintf('\n');
                end
            end
            
            function print_bar(nclasses)
                fprintf('      '); fprintf('%12s',repmat(' -----------',1,nclasses)); fprintf('\n');
            end
        end

        
        function cc_id        = identify_classes(MT)
            %
            % for AX = XB, if elements A and B are conjugate pairs for some other element X in the group,  they are in the same class
            %
            import am_lib.* am_dft.*

            nsyms = size(MT,2);
            % allocate space for conjugacy class
            cc_id = zeros(nsyms,1);
            % allocate space for conjugate elements
            conjugates = zeros(nsyms,1);
            % get inverse indicies
            I = (MT==1) * [1:size(MT,1)].';
            % determine conjugacy classes
            k = 0;
            for i = 1:nsyms
            if cc_id(i)==0
                k=k+1;
                % conjugate each element with all other group elements
                % A = X(j) * B * X(j)^-1
                for j = 1:nsyms
                    conjugates(j) = MT(j,MT(i,I(j)));
                end
                % for each subgroup element created by conjugation find the corresponding index of the element in the group
                % in order to save the class class_id number
                for j = 1:nsyms
                    cc_id( conjugates(j) ) = k;
                end
            end
            end
            % relabel classes based on how many elements each class has
            cc_id = reindex_using_occurances(cc_id);
        end

        function cart_vs_frac = get_symmetry_units(R)
            
            % if the inverse of every symmetry is its transpose, symmetries are in cartesian coordinates.
            for i = 1:size(R,3); iR(:,:,i) = inv(R(:,:,i)); end
            if all_(eq_(iR,permute(R,[2,1,3]))); cart_vs_frac = 'cart'; return; end
            
            % if there are just 0 and +/- 1s, it is definately in fractional coordinates [unless, of course, bas = eye(3)]
            if all(or(abs(R(:))==1,R(:)==0)); cart_vs_frac = 'frac'; return; end
            
            % otherwise
            error('unable to decipher units of symmetry');
        end

        function [bas]        = bas_hex2rhomb(bas)
            abc=am_cell.bas2abc(bas); c=abc(3); a=abc(1);
            abc(4:6) = acosd((2*c^2-3*a^2)/(2*c^2+6*a^2));
            abc(1:3) = sqrt(c^2+3*a^2)/3;
            
            a=abc(6);
            bas=abc(1)*[...
             cosd(a/2),-sind(a/2),0; 
             cosd(a/2),+sind(a/2),0;
             cosd(a)/cosd(a/2),0,sqrt(1-cosd(a)^2/cosd(a/2)^2)].';
        end
        
        function [bjc,gen]    = find_isomorphic_bijection(G,H)
            % bjc = find_isomorphic_bijection(g,h), g and h are multiplication tables 
            % isomorphism: g <==> h(:,:,inds) , returns inds = [] if groups are not isomorphic
            % % check with:
            % G = get_multiplication_table(g);
            % H = get_multiplication_table(h);
            % i=1; G - relabel_multiplication_table(H,bjc(:,i))

            import am_lib.* am_dft.*

            % number of symmetries
            nsyms = size(G,1);
            
            % sorting if desired (changes order of bjc but does not affect final result)
            Gfwd = [1:nsyms]; Hfwd = [1:nsyms];

            % get generators
            [gg,Gg] = get_generators(G); ngenerators = numel(gg);
            [hg,Hg] = get_generators_list(H,ngenerators);
            
            % compare multiplication tables, looking for bijections between symmetries g and h, i.e. isomorphism between G and H
            j=0; bjc=[]; gen=[];
            for i = 1:size(Hg,2)
                % fwd = relabeling of H to match G
                fwd(Gg) = Hg(:,i);
                if all(all( G == relabel_multiplication_table(H,fwd) ))
                    % isomorphism found,  save bijection which takes symmetries in H to G in original order
                    j=j+1; bjc(Gfwd,j) = Hfwd(fwd); gen(:,j) = hg(:,i);
                    % save
                end
            end
            nbjcs = j;

            % relabel generators for h in terms of g
            for j = 1:nbjcs
                bak(bjc(:,j)) = [1:nsyms]; gen(:,j) = bak(gen(:,j));
            end

        end

        function [fwd,M]      = find_pointgroup_transformation(g,h,algo)
            % finds the permutation fwd and rotation M which transforms point group g into h:
            % matmul_(matmul_(M(:,:,k),g(:,:,:)),inv(M(:,:,k))) == h(:,:,bjc(:,k)) for any k

            import am_dft.* am_lib.*
            
            % algo = 1, returns only 1 bijection of many possibilities
            if nargin<3; algo=0; end

            % find bijection bjc: h -> g and generators for h
            G=get_multiplication_table(g);
            H=get_multiplication_table(h);
            [bjc,gen] = find_isomorphic_bijection(G,H);
            srt_ = rankc_(bjc); bjc = bjc(:,srt_); gen = gen(:,srt_);

            % check diagonalize symmetry values to make sure the symmetries are of the correct type
            % this removes the possibility that multiplication table may be obeyed with symmetries of differnet types
            nbjcs=size(bjc,2); ex_ = true(1,nbjcs); 
            g_ps_id = identify_point_symmetries(g); 
            h_ps_id = identify_point_symmetries(h);
            for j = 1:nbjcs;   ex_ = all(h_ps_id(bjc(:,j))==g_ps_id,1); end
            bjc = bjc(:,ex_); gen = gen(:,ex_); nbjcs=size(bjc,2);

            % find integer transformation matrix (try all possible bijections and
            % integer transformations with nonzero-determinants); in reality, not all
            % matrices are tried, instead: 
            %       first try the identity
            %       then try all matrices with elements in [-1,0,1] with determinant 1 
            %       then try all matrices with elements in [-2,-1,0,1,2] with determinant 1 
            q = 0; ngens = size(gen,1); 
            for m = 0:1
                % generate the possible transformation matrices
                switch m
                    case 0
                        % first try the identity
                        X = eye(3); nXs = 1;
                    otherwise
                        % then try other integer matrices with nonzero determinants
                        N=9; Q=[-m:m]; nQs=numel(Q);[Y{N:-1:1}]=ndgrid(1:nQs); X=reshape(Q(reshape(cat(N+1,Y{:}),[],N)).',3,3,[]);
                        ex_ = false(1,size(X,3)); for i = 1:size(X,3); ex_(i) = abs(det(X(:,:,i)))<am_lib.eps; end
                        X = X(:,:,~ex_); nXs = size(X,3);
                        % sort it to make it nice (put positive values first)
                        X=X(:,:,rankc_( -[max(reshape(X,9,[]));min(reshape(X,9,[]));sum(reshape(X,9,[]),1)] ));
                end
                % get inverse elements
                invX = zeros(3,3,nXs); for i = 1:nXs; invX(:,:,i) = inv(X(:,:,i)); end

                % find similarity transform which converts all symmetries from g to h by
                % checking each matrix generated above
                for k = 1:nXs; for j = 1:nbjcs
                    gencheck_ = true;
                    % check generators
                    for i = 1:ngens
                        if any(any(abs( X(:,:,k) * g(:,:,gen(i,j)) * inv(X(:,:,k)) - h(:,:,bjc(gen(i,j),j)) )>am_lib.eps ))
                            gencheck_ = false; break;
                        end
                    end
                    % check whether a match has been found
                    if gencheck_
                        % at this point a bijection permutation and transformation matrix have been
                        % found such that: matmul_(matmul_(X(:,:,k),g(:,:,:)),invX(:,:,k)) == h(:,:,bjc(:,j))
                        % save it 
                        q = q+1;
                        M(:,:,q) = X(:,:,k); 
                        fwd(:,q) = bjc(:,j);
                        % only one is needed
                        if algo==1
                            return; 
                        end
                    end
                end;end
            end
        end
        
        
        % clusters

        % get_model('toy tb' , pc, cutoff, spdf)
        % get_model('toy bvk', pc, cutoff)
        % get_model('toy bvt', pc, cutoff)
        % get_model('tb' , pc, cutoff, spdf, dft, nskips)
        % get_model('bvk', pc, cutoff, uc, md)
        function [ip]         = get_model(flag,varargin)

            import am_lib.* am_dft.*
            
            % parse input
            if     contains(flag,'toy')
                if     contains(flag,'tb');  [pc, cutoff, spdf] = deal(varargin{:}); 
                elseif contains(flag,'bvk'); [pc, cutoff] = deal(varargin{:});
                end                
            else
                if     contains(flag,'tb');  [pc, cutoff, spdf, dft, nskips] = deal(varargin{:}); 
                elseif contains(flag,'bvk'); [pc, cutoff, uc, md] = deal(varargin{:});
                end
            end
            
            % set variable names
            if     contains(flag,'tb');  M_ = 'v_tb'; ME_='sk'; N_ = 'nbands';    H_F_ = 'H'; H_C_ = 'Hc';
            elseif contains(flag,'bvk'); M_ = 'v_fc'; ME_='fc'; N_ = 'nbranches'; H_F_ = 'D'; H_C_ = 'Dc'; end

            % get irreducible pairs
            fprintf(' ... identifying irreducible pairs'); tic;
            ip = get_irreducible_cluster(pc,2,cutoff);
            fprintf(' (%.f secs)\n',toc);
            
            % record model
            ip.model = flag;

            % print irreducible pair cluster info 
            print_cluster(ip,'bond,cart');
            
            % get irreducible shells
            fprintf(' ... analyzing matrix element symmetries'); tic;
            if     contains(flag,'tb');  [ip.W,ip.(M_),ip.Dj] = get_cluster_matrix_elements(ip,spdf);
            elseif contains(flag,'bvk'); [ip.W,ip.(M_),ip.Dj] = get_cluster_matrix_elements(ip); end
            fprintf(' (%.f secs)\n',toc);

            % get hamiltonian
            fprintf(' ... building Hamiltonian'); tic;
            H = get_cluster_hamiltonian(ip); ip.(N_) = size(H,1); ip.(H_C_) = H; ip.(H_F_) = matlabFunction(ssum_(ip.(H_C_))); 
            fprintf(' (%.f secs)\n',toc);

            if ~contains(flag,'toy')
                % tight binding model
                fprintf(' ... solving for matrix element values '); tic;
                if     contains(flag,'tb');  ip.(ME_) =  get_tb_parameters(ip,dft,nskips);
                elseif contains(flag,'bvk'); ip.(ME_) = get_bvk_parameters(ip,uc,md); 
                end
                fprintf(' (%.f secs)\n',toc);

%                 if  contains(flag,'bvk')
%                     fprintf(' ... enforcing acoustic sum rules'); tic;
%                     fc_before = [ip.fc{:}]; [ip.fc] = enforce_asr(ip); fc_after = [ip.fc{:}]; 
%                     fprintf(', max difference = %.2f%%',max(abs(fc_before-fc_after)./abs(fc_before))*100);
%                     fprintf(' (%.f secs)\n',toc);
%                 end

            end
            
            function [tb]          = get_tb_parameters(ip,dft,nskips)
                % nskips : number of dft bands to skip (e.g. 5)

                % fit neighbor parameter at high symmetry points using poor man's simulated anneal
                d = am_lib.normc_(ip.bas*(ip.tau(:,ip.cluster(2,:))-ip.tau(:,ip.cluster(1,:))));
                d4fc = repelem(d,cellfun(@(x)size(x,2),ip.W)); nfcs=numel(d4fc); x=zeros(1,nfcs);
                d=unique(am_lib.rnd_(d4fc)); d=conv([d,Inf],[1 1]/2,'valid'); nds = numel(d); r_best = Inf;

                % set simulated annealing temeprature and optimization options
                kT = 20; kT_decay_ = @(kT,i) kT .* exp(-i/5); rand_ = @(x) (0.5-rand(size(x))).*abs(x./max(x));
                opts = optimoptions('lsqnonlin','Display','None','MaxIter',7);

                % select bands
                bnd_id = [1:ip.nbands]+nskips;

                % define cost function
                kpt_id = 1:max(round(dft.nks/20),1):dft.nks;
                % kpt_id = [1:dft.nks];
                cost_ = @(x) dft.E(bnd_id,kpt_id) - am_dft.eval_energies_(ip,x,dft.k(:,kpt_id));

                % poor man's simulated annealing: loop over distances, incorporating each shell at a time
                % it appears that ignoring the loop over distance is better, at least for cases with small pair cutoffs
                for j = 2:nds
                    for i = 1:30
                        % simulated annealing
                        if i ~= 1; x = x_best + rand_(x_best) * kT_decay_(kT,i); end

                        % optimize
                        [x,r] = am_lib.lsqnonlin_(cost_, x, [d4fc>d(j)], [], [], opts);

                        % save r_best parameter
                        if r < r_best; r_best = r; x_best = x;
                            % plot band structure (quick and dirty)
                            plot([1:dft.nks], am_dft.eval_energies_(ip,x,dft.k),'-k',...
                                 [1:dft.nks], dft.E(bnd_id,:),':r');
                            set(gca,'XTick',[]); axis tight; grid on;
                            ylabel('Energy E'); xlabel('Wavevector k'); drawnow;
                        end
                    end
                end

                % redefine cost function on all kpoints
                kpt_id = [1:dft.nks];
                cost_ = @(x) dft.E(bnd_id,kpt_id) - am_dft.eval_energies_(ip,x,dft.k(:,kpt_id));

                % final pass with all parameters and all kpoints
                [x,~] = am_lib.lsqnonlin_(cost_,x,false(1,nfcs),[],[],opts);

                % save refined matrix elements and conform to ip
                for i = [1:ip.nclusters]; d(i)=size(ip.W{i},2); end; Evsk=cumsum(d); Svsk=Evsk-d+1;
                for i = [1:ip.nclusters]; tb{i} = x(Svsk(i):Evsk(i)); end
            end

            function [bvk]         = get_bvk_parameters(ip,uc,md)
                % Extracts symmetry adapted force constants.
                %
                % Q: What is the difference of these methods?
                % A: Force constants extracted, dispersions look different:
                %    Method 1:    8.9932    0.2321   -0.4806   -0.1065   -0.0337   -4.9099   14.6953   -0.5217    0.6199   -0.2590
                %    Method 2:   10.0033    0.4073   -0.5706    0.0775   -0.0314   -5.3528   16.0851   -0.5830    0.2604   -0.4623
                % Q: How do the methods comapre in terms of memory, speed, and
                %    ability to reproduce AIMD forces?
                % A:        METHOD     1           2
                %              R^2     0.668       0.633
                %      MATRIX SIZE     3*m*nFCs    3*n*m
                %          TIME [s]    0.322       0.031
                %      TIMES SLOWER    10x         1x

                import am_lib.* am_dft.*

                % [cart] get displacements and forces
                u = am_lib.matmul_( md.bas, am_lib.mod_( md.tau-uc.tau +.5 )-.5 );
                f = am_lib.matmul_( md.bas, md.force );
                
                % covert symmetries [pc-frac] to [cart] -- very important!
                sym_rebase_ = @(B,S) [[ am_lib.matmul_(am_lib.matmul_(B,S(1:3,1:3,:)),inv(B)), ...
                    reshape(am_lib.matmul_(B,S(1:3,4,:)),3,[],size(S,3))]; S(4,1:4,:)];
                Q_cart = ip.Q; Q_cart{1} = sym_rebase_(ip.bas, Q_cart{1}); Q_cart{1} = am_lib.wdv_(Q_cart{1});

                % select a method (method 1 is the most accurate)
                switch 2
                    case 1
                        % Method incorporating full intrinsic and crystalographic symmetries
                        % Note that ASR are automatically enforced when this method is used!
                        % F [3m * 1] = - U [3m * nFcs] * FC [ nFCs * 1]: n pairs, m atoms ==> FC = - U \ F
                        %
                        % Q: How are roto-inversions incorporated?
                        % A: Check with:
                        %
                        %         % define input prep
                        %         prep_ = @(x) [x(1),x(2),x(3),0,0,0,0,0,0;0,0,0,x(1),x(2),x(3),0,0,0;0,0,0,0,0,0,x(1),x(2),x(3)];
                        %         % select irreducible pair j and symmetry i
                        %         j=2;i=9;
                        %         % define symbolic displacement u and force constants c
                        %         u=sym('u_%d',[3,1]); c=sym('c_%d',[size(ip.W{j},2),1]);
                        %         % get rotations and inverse
                        %         R=pp.Q{1}(1:3,1:3,pp.q{1}(i)); iR=pp.Q{1}(1:3,1:3,pp.iq{1}(i));
                        %         % check without rotation
                        %         prep_(u)*ip.W{j} - equationsToMatrix(reshape(ip.W{j}*c,3,3)*u,c)
                        %         % check with rotation
                        %         R*prep_(iR*u)*ip.W{j} - equationsToMatrix(R*reshape(ip.W{j}*c,3,3)*iR*u,c)
                        %
                        % get U matrix and indexing I
                        [U,I] = am_dft.get_U_matrix(ip, uc, u);

                        % solve for force constants
                        fc = - U \ f(I);

                        % get number of force constants per shell
                        s_id = repelem([1:ip.nclusters],cellfun(@(x)size(x,2),ip.W));
                        
                        % save force constants
                        for s = 1:ip.nclusters; bvk{s} = double(fc(s_id==s).'); end
                    case 2
                        % basic method using full 3x3 second-order tensor (ignores intrinsic force constant symmetry)
                        % F [3 * m] = - FC [3 * 3n] * U [3n * m]: n pairs, m atoms ==> FC = - F / U
                        pp = am_dft.get_primitive_cluster(ip,uc); 

                        % construct force constant matrices
                        phi = []; pc_natoms = numel(pp.c_id);
                        for m = 1:pc_natoms 
                            % get forces : f = [ (x,y,z), (1:natoms)*nsteps ]
                            % get displacements : u [ (x,y,z)*orbits, (1:natoms)*nsteps ]
                            %  ... and solve for the generalized force constants: FC = - f / u
                            fc = - reshape( f(:,pp.c_id{m},:) , 3                , pp.ncenters(m) * md.nsteps) /...
                                   reshape( u(:,pp.o_id{m},:) , 3 * pp.norbits(m), pp.ncenters(m) * md.nsteps);
                            % reshape into 3x3 matrices
                            phi = double(cat(3,phi,reshape(fc,3,3,[])));
                        end

                        % reorder force constants to match pp (phi is ordered in terms of m 
                        % here, but pp is in terms of irreducible pairs)
                        phi(:,:,cat(1,pp.pp_id{:})) = phi;

                        % transform force constants according to symmetries Q (used to go between orbit and irrep)
                        q = cat(1,pp.q_id{:});
                        for j = 1:size(phi,3); phi(:,:,j) = permute( phi(:,:,j), Q_cart{2}(:,q(j)) ); end
                        phi = am_lib.matmul_( am_lib.matmul_( Q_cart{1}(1:3,1:3,q) ,phi), permute(Q_cart{1}(1:3,1:3,q),[2,1,3]) );

                        % solve for symmetry adapted force constants : A x = B
                        for i = 1:ip.nclusters
                            % find primitive clusters involving irreducible cluster i
                            ex_ = pp.pp2ip==i;
                            if any(ex_)
                                A = repmat(double(ip.W{i}),sum(ex_),1);
                                B = reshape(phi(:,:,ex_),[],1);
                                % get force constants as row vectors
                                bvk{i} = reshape( A \ B , 1, []);
                            end
                        end
                end
            end
            
            function [W,M,Dj]      = get_cluster_matrix_elements(ip,spdf)
                % [ip ] = analyze_cluster_matrix_elements(flag,ip,spdf)
                %
                % set orbitals per irreducible atom: spdf = {'d','p'};
                % may wish to do this to set it per species: x={'p','d'}; spdf={x{ic.species}};
                %
                % Q: What should zeroth, first, and second neighbor force constants look
                %    like for a typical material?
                % A: The force constant matrices for Si are listed below:
                %
                %     zeroth neighbor [0.00000    0.00000    0.00000]
                % 
                %     [ c01_11,      0,      0]
                %     [      0, c01_11,      0]
                %     [      0,      0, c01_11]
                % 
                %     first neighbor  [1.36750    1.36750    1.36750]
                % 
                %     [ c03_11, c03_21, c03_21]
                %     [ c03_21, c03_11, c03_21]
                %     [ c03_21, c03_21, c03_11]
                % 
                %     second neighbor [2.73500    2.73500    0.00000]
                % 
                %            (INCORRECT)	         (CORRECT, MELVIN LAX p 264)
                %     [ c02_11, c02_21,      0]  ==>  [ c02_11, c02_21, c02_31]
                %     [ c02_21, c02_11,      0]  ==>  [ c02_21, c02_11, c02_31]
                %     [      0,      0, c02_33]  ==>  [-c02_31,-c02_31, c02_33]
                %
                %     The erroneous version of this force constant matrix was published in 
                %
                %       H. M. J. Smith, Philosophical Transactions of the Royal Society a:
                %       Mathematical, Physical and Engineering Sciences 241, 105 (1948).  
                %
                %     as pointed out by
                %
                %   	M. Lax, Symmetry Principles in Solid State and Molecular Physics
                %   	(Dover Publications, 2012), p 264 and p 364.
                %
                %     Another publication with the corrected version of the force constant
                %     is:
                %
                %       R. Tubino and J. L. Birman, Phys. Rev. B 15, 5843 (1977).
                %
                %     Interestingly, that the second-neighbor force constant is
                %     anti-symmetric implies that the order of differentiation matters and
                %     the intrinsic transpositional symmetry does not apply. 
                %
                %     And, lastly, ...
                %
                %     third neighbor  [4.06600    1.35500    1.35500]
                % 
                %     [ c04_11, c04_21, c04_21]
                %     [ c04_21, c04_22, c04_32]
                %     [ c04_21, c04_32, c04_22]
                %
                %
                % Q: How are symmetry operations represented in the flattened basis?
                % A: Check it out with the code below.
                %
                %     0) initialize symbolic variables
                %        a=sym('a_%d%d',[3,3]);
                %        b=sym('b_%d%d',[3,3]);
                %        c=sym('c_%d%d',[3,3]);
                % 
                %     1) matrix multiplication 
                %        equationsToMatrix(a*b*c.',b(:)) - kron(c,a)
                % 
                %     2) flip symmetry
                %        xx = [1:9]; xt = reshape(permute(reshape(xx,3,3),ip.Q{2}(:,iq)),1,9);
                %        T = zeros(9,9); T(sub2ind([9,9],xx,xt) = 1;
                %        equationsToMatrix(a*(b.')*c.',b(:)) - kron(c,a) * T
                % 
                %     3) flip symmetry (another way)
                %        X=equationsToMatrix(a*b*c.',b(:)); 
                %        simplify(equationsToMatrix(a*(b.')*c.',b(:)) - X(:,[1,4,7,2,5,8,3,6,9]))
                % 
                % --------------------------------------------------------------------------------
                %
                % Q: What is the intrinsict symmetry of third-order force
                %    constants due to permutation of coordinate axes x,y,z?
                % A: Each page of the tensor is symmetry.
                %
                %     phi(:,:,1) =				   phi(1,:,:) =
                %         [ c_111, c_211, c_311]   		[ c_111, c_211, c_311]
                %         [ c_211, c_221, c_321]   		[ c_211, c_221, c_321]
                %         [ c_311, c_321, c_331]   		[ c_311, c_321, c_331]
                %     phi(:,:,2) =				   phi(2,:,:) =
                %         [ c_211, c_221, c_321]   		[ c_211, c_221, c_321]
                %         [ c_221, c_222, c_322]   		[ c_221, c_222, c_322]
                %         [ c_321, c_322, c_332]   		[ c_321, c_322, c_332]
                %     phi(:,:,3) =				   phi(3,:,:) =
                %         [ c_311, c_321, c_331]   		[ c_311, c_321, c_331]
                %         [ c_321, c_322, c_332]   		[ c_321, c_322, c_332]
                %         [ c_331, c_332, c_333]   		[ c_331, c_332, c_333]
                %
                %    And the same thing for the middle phi(:,1:3,:). Check if out with the code:
                %
                %         % get [x,y,z] flip/permutation operators
                %         a = reshape([1:27],3,3,3); p = perms(1:3).'; nps = size(p,2); F = zeros(27,27,nps);
                %         for i = 1:nps; F(:,:,i) = sparse( flatten_(a), flatten_(permute(a,p(:,i))), ones(27,1), 27, 27); end
                %         % enforce intrinsic symmetry (immaterial order of differentiation: c == c.')
                %         W = sum(F-eye(27),3);
                %         % get linearly-independent nullspace and normalize to first nonzero element
                %         W=real(null(W)); W=frref_(W.').'; W(abs(W)<am_lib.eps)=0; W(abs(W-1)<am_lib.eps)=1; W=W./accessc_(W,findrow_(W.').');
                %         % define parameters
                %         c = sym(sprintf('c_%%d%%d%%d',i),[3,3,3],'real'); c = c(findrow_(double(W).'));
                %         % get symmetry adapted force constants
                %         phi = reshape( sym(W)*c(:), [3,3,3]);
                %
                % Q: How do second and third rank tensors rotate?
                % A:
                %         %% rotate 2nd rank tensor A
                %         clear;clc;rng(1); R=rand(3,3); A = rand(3,3);
                %         AR = zeros(3,3);
                %         for m = 1:3; for n = 1:3;
                %         for i = 1:3; for j = 1:3
                %             AR(m,n) = AR(m,n) + A(i,j)*R(m,i)*R(n,j);
                %         end; end
                %         end; end
                %         AR - R*A*R.'
                %         AR - reshape(kron(R,R)*A(:),3,3)
                %
                %         %% rotate 3rd rank tensor A
                %         clear;clc;rng(1); R=rand(3,3); A = rand(3,3,3);
                %         AR = zeros(3,3,3);
                %         for m = 1:3; for n = 1:3; for o = 1:3
                %         for i = 1:3; for j = 1:3; for k = 1:3
                %             AR(m,n,o) = AR(m,n,o) + A(i,j,k)*R(m,i)*R(n,j)*R(o,k);
                %         end; end; end
                %         end; end; end
                %         AR - reshape(kron(kron(R,R),R)*A(:),[3,3,3])
                %
                % Q: How do triple dot products work? For example, when
                %    fourier-transforming third-order force constants.
                % A: The triple dot product PHI ... u1 u2 u3, as described by Ziman,
                %    produces as a scalar and is given in Einstein summation as
                %    PHI(ijk) u1(i) u2(j) u3(k).
                %
                %         % apply double dot product
                %         F = zeros(3,1);
                %         for m = 1:3
                %         for i = 1:3; for j = 1:3
                %             F(m) = F(m) + AR(i,j,m)*u(i)*w(j);
                %         end; end
                %         end
                %
                % Q: What about a double dot product? For example, when
                %    third-order force constants multiply displacements to
                %    obtain forces.
                % A: Check it out with the code below:
                %
                %         % double dot product of rotated 3rd rank tensor
                %         clear;clc;rng(1); R=rand(3,3); u = rand(3,1); w = rand(3,1); A = rand(3,3,3);
                %         % rotate A
                %         AR = zeros(3,3,3);
                %         for m = 1:3; for n = 1:3; for o = 1:3
                %         for i = 1:3; for j = 1:3; for k = 1:3
                %             AR(m,n,o) = AR(m,n,o) + A(i,j,k)*R(m,i)*R(n,j)*R(o,k);
                %         end; end; end
                %         end; end; end
                %         % apply double dot product
                %         F = zeros(3,1);
                %         for m = 1:3
                %         for i = 1:3; for j = 1:3
                %             F(m) = F(m) + AR(i,j,m)*u(i)*w(j);
                %         end; end
                %         end
                %         % equivalent formulations incorproating rotation:
                %         F - reshape(kron(kron(R,R),R)*A(:),9,3).'*reshape(u(:)*w(:).',[],1)
                %         % equivalent formulations without incorporating rotation
                %         F - reshape(          AR          ,9,3).'*reshape(u(:)*w(:).',[],1)
                %         F - [ u.'*AR(:,:,1)*w;  u.'*AR(:,:,2)*w ; u.'*AR(:,:,3)*w ]
                %         F - matmul_(AR,w).'*u
                %         F - reshape(AR,9,3).'*flatten_(u*w.')
                %         F - squeeze(sum(sum( AR.*(u*w.') ,1),2))

                digits(10); 

                % if spdf is suppied, tight binding calculation, else force constant calculation
                if nargin<2; spdf=repmat({'p'},1,max(ip.x2i)); end

                % get stabilizers for each irreducible cluster
                [~,~,~,s_ck] = am_dft.get_action(ip);

                % reversal symmetries should not be ignored, however in the rare even that
                % one would want to ignore them, set ignore_reversal_symmerties to true.
                ignore_reversal_symmerties = false;
                if ignore_reversal_symmerties; s_ck((end/2):end,:)=false; end

                % covert symmetries [pc-frac] to [cart] -- very important!
                sym_rebase_ = @(B,S) [[ am_lib.matmul_(am_lib.matmul_(B,S(1:3,1:3,:)),inv(B)), ...
                    reshape(am_lib.matmul_(B,S(1:3,4,:)),3,[],size(S,3))]; S(4,1:4,:)];
                Q_cart = ip.Q; Q_cart{1} = sym_rebase_(ip.bas, Q_cart{1});

                % for each irreducible atom, set azimuthal quantum numbers J{:}, symmetries D{:}, and parity-transpose F{:}
                % for force constants, vectors transform like 'p' orbitals and all J = 3.
                Dj = get_tb_symmetry_representation(spdf, Q_cart);

                % get dimensions
                d = cellfun(@(x)size(x,1),Dj);
                % get form of force constants for irreducible prototypical bonds
                for p = 1:ip.nclusters
                    % get irreducible atom indicies
                    ijk = ip.x2i( ip.cluster(:,p) ); 
                    % use stabilzer group to determine crystallographic symmetry relations; A*B*C' equals kron(C,A)*B(:)
                    % a=sym('a',[3,3]); b=sym('b',[3,3]); c=sym('b',[3,3]); simplify(a*b*c.'-reshape(kron(c,a)*b(:),3,3))
                    % a=sym('a',[3,3]); b=sym('b',[3,5]); c=sym('c',[5,5]); simplify(a*b*c.'-reshape(kron(c,a)*b(:),3,5))
                    % choose an algorithm ALGO=2 is the default and correct version
                    switch 2
                        case 1
                            % original (and erroneous version) just like in Smith's paper on force constants of diamond.
                            W{p} = ones(1,1,sum(s_ck(:,p)));
                            for j = ijk
                                W{p} = am_lib.kron_( Dj{j}(:,:,s_ck(:,p)) , W{p} );
                            end
                        case 2
                            % incoprorating flip symmetry when atomic indicies are flipped, the matrix element is transposed. 
                            % T=zeros(9,9); T(sub2ind([9,9],[1:9],reshape(reshape([1:9],3,3).',1,9)))=1;
                            % a=sym('a',[3,3]); b=sym('b',[3,5]); c=sym('c',[5,5]); equationsToMatrix(a*b*c.',b(:))     - kron(c,a)
                            % a=sym('a',[3,3]); b=sym('b',[3,3]); c=sym('c',[3,3]); equationsToMatrix(a*(b.')*c.',b(:)) - kron(c,a) * T  % FIRST TRANSPOSE THEN MULTIPLY
                            % a=sym('a',[3,3]); b=sym('b',[3,3]); c=sym('c',[3,3]); equationsToMatrix( (a*b*c.').',b(:)) - T * kron(c,a) % FIRST MULTIPLY THEN TRANSPOSE
                            sym_list = find(s_ck(:,p)); W{p} = ones(1,1,sum(s_ck(:,p)));
                            for j  = ijk;    W{p} = am_lib.kron_( Dj{j}(:,:,sym_list) , W{p} ); end
                            for wi = 1:numel(sym_list); W{p}(:,:,wi) = am_lib.transpose_( d(ijk), Q_cart{2}(:,sym_list(wi)) ) * W{p}(:,:,wi); end
                    end
                    % make equations
                    W{p} = sum( W{p} - eye(prod(d(ijk))), 3);
                    % get linearly-independent nullspace and normalize to first nonzero element
                    W{p} = null(W{p}); W{p} = am_lib.frref_(W{p}.').'; W{p} = W{p}./am_lib.accessc_(W{p},am_lib.findrow_(W{p}.').'); W{p} = am_lib.wdv_(W{p}); 
                    % define parameters
                    C{p} = sym( sprintf('c%02i%s',p,repmat('_%d',1,ip.nvertices)), d(ijk), 'real'); C{p} = C{p}(am_lib.findrow_(W{p}.'));
                    % this must come before the sort below!
                    M{p} = reshape( sym(W{p})*C{p}(:), d(ijk) );
                    % reorder C to be in the same order as D
                    [C{p},n] = sort(C{p}(:).'); W{p} = W{p}(:,n); 
                end

                function Dj = get_tb_symmetry_representation(spdf,Q_cart)
                    % set symmetries D{:} per irreducible atom
                    % irreducible atom given a list of orbitals for each
                    % irreducible atom, spdf = {'sp','d'}

                    % get symmetries
                    nQs = size(Q_cart{2},2);

                    % transform symmetries to the tight binding representation (wiger functions)
                    W_=cell(1,3); for j_=[1:3]; W_{j_} = am_lib.get_wigner(j_,Q_cart{1}(1:3,1:3,:),'real'); end

                    % set orbitals J{:}, symmetries D{:}, and parity-transpose T{:} for each irreducible atom
                    ic_natoms = numel(spdf); Dj = cell(1,ic_natoms);
                    for i_ = 1:ic_natoms
                        % set orbitals
                        J{i_} = am_lib.findrow_('spdf'==spdf{i_}(:)).'-1;
                        % set start and end points for J
                        E=cumsum(J{i_}*2+1); S=E-(J{i_}*2+1)+1;
                        % construct D matrix
                        d_ = max(E); Dj{i_} = zeros(d_,d_,nQs);
                        for j_ = 1:length(J{i_})
                            if J{i_}(j_)==0 % s orbital
                                Dj{i_}(S(j_):E(j_),S(j_):E(j_),:) = 1;
                            else % p,d,f orbitals
                                Dj{i_}(S(j_):E(j_),S(j_):E(j_),:) = W_{J{i_}(j_)};
                            end
                        end
                        % correct rounding error
                        Dj{i_} = am_lib.wdv_(Dj{i_});
                    end
                end
            end
            
            function [H]           = get_cluster_hamiltonian(ip)
                %
                % Q: What do the force constant matrices loop like when rotated from the
                %    irreducible bond to the orbit? 
                % A: Clusters up to first-neighbor in Si, indexed by the pair of atoms are:
                %
                %    pp.cluster = 
                %      1     2     1     1     1     1     2     2     2     2
                %      1     2     2     3     6     7     1     4     5     8
                %
                %    The force constant matrices for each first-neighbor orbit involving
                %    the primitive cell atom at [0 0 0] are:
                %
                %     for p = 3									for p = 4
                %             rij = 										rij = 
                %                   [ 1.3552, 1.3552, 1.3552]					  [ 1.3552,-1.3552,-1.3552]
                %             matrix_ = 									matrix_ = 
                %                   [ c02_11, c02_21, c02_21]					  [ c02_11,-c02_21,-c02_21]
                %                   [ c02_21, c02_11, c02_21]					  [-c02_21, c02_11, c02_21]
                %                   [ c02_21, c02_21, c02_11]					  [-c02_21, c02_21, c02_11]
                %
                %     for p = 5									for p = 6
                %             rij = 										rij = 
                %                   [-1.3552, 1.3552,-1.3552]					  [-1.3552,-1.3552, 1.3552]
                %             matrix_ = 									matrix_ = 
                %                   [ c02_11,-c02_21, c02_21]					  [ c02_11, c02_21,-c02_21]
                %                   [-c02_21, c02_11,-c02_21]					  [ c02_21, c02_11,-c02_21]
                %                   [ c02_21,-c02_21, c02_11]					  [-c02_21,-c02_21, c02_11]
                %
                %    For the primitive cell atom at [1 1 1]/4, the force constant matrices are:
                %   
                %     for p = 7									for p = 10
                %             rij = 										rij = 
                %                   [-1.3552,-1.3552,-1.3552]					  [-1.3552, 1.3552, 1.3552]
                %             matrix_ = 									matrix_ = 
                %                   [ c02_11, c02_21, c02_21]					  [ c02_11,-c02_21,-c02_21]
                %                   [ c02_21, c02_11, c02_21]					  [-c02_21, c02_11, c02_21]
                %                   [ c02_21, c02_21, c02_11]					  [-c02_21, c02_21, c02_11]
                %
                %     for p = 9 								for p = 8
                %             rij = 										rij = 
                %                   [ 1.3552,-1.3552, 1.3552]					  [ 1.3552, 1.3552,-1.3552]
                %             matrix_ = 									matrix_ = 
                %                   [ c02_11,-c02_21, c02_21]					  [ c02_11, c02_21,-c02_21]
                %                   [-c02_21, c02_11,-c02_21]					  [ c02_21, c02_11,-c02_21]
                %                   [ c02_21,-c02_21, c02_11]					  [-c02_21,-c02_21, c02_11]
                %
                %   Compare these to Melvin Lax p 266.

                digits(10); 

                % switch between phonon and electrons
                if     contains(ip.model,'tb');  matrix_element_ = 'v_tb'; 
                elseif contains(ip.model,'bvk'); matrix_element_ = 'v_fc'; end

                % expand cluster
                pp = am_dft.get_primitive_cluster(ip);

                % get hamiltonian block dimensions and start/end sections
                d(ip.x2p(ip.cluster(1,:))) = cellfun(@(x)size(x,1),ip.(matrix_element_)); 
                d(ip.x2p(ip.cluster(2,:))) = cellfun(@(x)size(x,2),ip.(matrix_element_));
                E=cumsum(d); S=E-d+1; nbands=E(end);

                % construct symbolic Hamiltonian matrix
                H = sym(zeros(nbands,nbands,ip.nclusters)); kvec = sym('k%d',[3,1],'real'); 
                for p = 1:pp.nclusters
                    %    c     = irreducible cluster index
                    %    m,n   =   primitive cell atomic indicies
                    c = pp.pp2ip(p); 
                    m = pp.x2p(pp.cluster(1,p)); mp = S(m):E(m);
                    n = pp.x2p(pp.cluster(2,p)); np = S(n):E(n); 

                    % get bond vector
                    rij = pp.tau(:,pp.cluster(2,p)) - pp.tau(:,pp.cluster(1,p)); rij = am_lib.wdv_(rij,am_dft.tiny);

                    % fourier transform build hamiltonian matrix
                    H(mp,np,c) = H(mp,np,c) + pp.(matrix_element_){p} .* exp(sym(2i*pi) * sym(rij(:).') * kvec(:) );
                end

                % simplify (speeds evaluation up significantly later)
                for i = 1:ip.nclusters; H(:,:,i) = simplify(rewrite(H(:,:,i),'cos'),'steps',3); end

                % multiply mass factor to get the dynamical matrix
                if contains(ip.model,'bvk')
                    % get unique masses
                    mass = sym('m_%d',[1,numel(unique(ip.species))],'real');
                    % get mass matrix
                    p2i=[]; p2i(ip.x2p)=ip.x2i; mass=mass(repelem(p2i,1,3)); mass=(mass.'*mass); 
                    % include mass
                    for i = 1:ip.nclusters; H(:,:,i) = H(:,:,i)./sqrt(mass); end
                end
            end
            
            function                 print_cluster(ip,flags)

                if nargin<2; flags=''; end

                if isfield(ip,'p2i'); isprimitive=true;  str_type='primitive';
                else;                 isprimitive=false; str_type='irreducible'; end

                % choose between fractional and cartesian coordinates
                if contains(flags,'frac') || ~contains(flags,'cart'); ip.bas=eye(3); end

                fprintf(' ... %s %i-atom clusters:\n', str_type, ip.nvertices);

                fprintf('     '); 
                    fprintf(['%-',num2str(4+3*ip.nvertices),'s'],'cluster');
                    if contains(flags,'bond')
                        fprintf('%24s',sprintf('atom #%i',1)); for i = 2:ip.nvertices; fprintf('%24s',sprintf('bond (1-%i)',i)); end
                    else
                        for i = 1:ip.nvertices; fprintf('%24s',sprintf('atom #%i',i)); end
                    end
                    fprintf('%8s','length'); 
                    if isprimitive; fprintf('%5s','p2i'); end
                fprintf('\n');  

                % bar
                fprintf('     '); 
                    fprintf(['%-',num2str(4+3*ip.nvertices),'s'],repmat('-',4+3*ip.nvertices,1));
                    for i = 1:ip.nvertices; fprintf('%8s%8s%8s','-------','--------','--------'); end
                    fprintf('%8s','-------'); 
                    if isprimitive; fprintf('%5s','----'); end
                fprintf('\n');  

                % table entries
                switch ip.nvertices
                    case 2
                        REF = ip.bas*ip.tau(:,ip.cluster(1,:));
                        POS = ip.bas*ip.tau(:,ip.cluster(2,:));
                        d   = am_lib.normc_(REF-POS);
                    case 3
                        REF = reshape(ip.bas*ip.tau(:,ip.cluster(1,:)),[3,1           ,ip.nclusters]);
                        POS = reshape(ip.bas*ip.tau(:,ip.cluster     ),[3,ip.nvertices,ip.nclusters]);
                        d   = am_lib.normc_(REF-POS)/sqrt(ip.nvertices);
                    otherwise 
                        error('not yet implemented');
                end

                for i = 1:ip.nclusters
                fprintf('     '); 
                    % cluster number and symbols
                    fprintf('%-4i %2s',i,ip.symb{ip.species(ip.cluster(1,i))}); fprintf('-%2s',ip.symb{ip.species(ip.cluster(2:end,i))}); 
                    % center and bonds
                    if contains(flags,'bond')
                        fprintf(' %7.3f', ip.bas*ip.tau(:,ip.cluster(1,i)) ); 
                        fprintf(' %7.3f', ip.bas*ip.tau(:,ip.cluster(2:end,i))-ip.bas*ip.tau(:,ip.cluster(1,i))); 
                    else
                        fprintf(' %7.3f', ip.bas*ip.tau(:,ip.cluster(:,i)));
                    end
                    % length
                    fprintf(' %7.3f',d(i)); 
                    % irreducible index if primitive
                    if isprimitive; fprintf('%5i', ip.pp2ip(i)); end
                fprintf('\n');
                end
                
                % print orbit and stabilizer group generators
                if ~isprimitive
                    % get stabilizers and generators in fractional coordinates before performing change of basis
                    [~,~,~,s_ck,g_ck] = am_dft.get_action(ip);
                    
                    % print only generators of stabilizers and orbits
                    Q_ex_ = @(ex_) deal(ip.Q{1}(:,:,ex_),ip.Q{2}(:,ex_));
                    for i = 1:ip.nclusters
                        [Q_{1:2}] = Q_ex_(s_ck(:,i)); tmp = am_dft.get_generators(am_dft.get_multiplication_table(Q_)); s_ck(:,i) = false; s_ck(tmp,i) = true;
                    end
                    
                    % convert symmetries if necessary
                    sym_rebase_ = @(B,S) [[ am_lib.matmul_(am_lib.matmul_(B,S(1:3,1:3,:)),inv(B)), ...
                                    reshape(am_lib.matmul_(B,S(1:3,4,:)),3,[],size(S,3))]; S(4,1:4,:)];

                    % get field properly
                    if     isfield(ip,'S'); sym='S'; ip.(sym)    = sym_rebase_(ip.bas,ip.(sym));    ip.(sym)    = am_lib.wdv_(ip.(sym));
                    elseif isfield(ip,'Q'); sym='Q'; ip.(sym){1} = sym_rebase_(ip.bas,ip.(sym){1}); ip.(sym){1} = am_lib.wdv_(ip.(sym){1}); end
                    
                    % get names
                    sym_name = am_dft.get_long_ss_name(ip.(sym));

                    % print generators and stabilizers (bond group + reversal group)
                    for i = 1:ip.nclusters%; if am_lib.gt_(d(i),0)
                        fprintf('     ----------------------------------------------------- cluster #%-3i\n', i);
                        print_helper_(sym_name,'generators',g_ck(:,i));
                        try
                        print_helper_(sym_name,'stabilizers',s_ck(:,i));
                        catch
                            asdf
                        end
                        fprintf('\n');
                    end%;end
                end

                function print_helper_(sym_name,name,x_)
                    nGs = sum(x_); nentries_ = 2; print_length = 2+max(cellfun(@(x)numel(strtrim(x)),sym_name));
                    if nGs>0
                    fprintf('     %-12s',[name,':']);
                        for j = 1:nGs
                            G=sym_name(x_);
                            fprintf('%*s',print_length,strtrim(G{j})); 
                            if mod(j,nentries_)==0 && j~=numel(G)
                                fprintf('\n'); 
                                fprintf('     %-10s  ',' '); 
                            end
                        end
                        fprintf('\n');
                    end
                end
            end

            function [bvk]         = enforce_asr(ip)

                % expand parameters over primitive clusters as well
                pp = am_dft.get_primitive_cluster(ip);

                % get self-force constants
                uc_natoms = max(pp.x2p);
                phi = zeros(3,3,uc_natoms);
                for m = 1:uc_natoms
                for i = 1:pp.nclusters
                    if m == pp.x2p(pp.cluster(1,i))   &&   pp.cluster(1,i) ~= pp.cluster(2,i)
                        phi(:,:,m) = phi(:,:,m) + pp.fc{i};
                    end
                end
                end
                phi = - phi;

                % copy 
                bvk = ip.fc;
                
                % solve for symmetry-adapted self-force constants : A x = B
                for m = 1:uc_natoms
                for i  = 1:ip.nclusters
                    if m == ip.x2p(ip.cluster(1,i))   &&   ip.cluster(1,i) == ip.cluster(2,i)
                        A = double(ip.W{i});
                        B = reshape(phi(:,:,m),[],1);
                        % get force constants as row vectors
                        bvk{i} = reshape( A \ B , 1, []);
                    end
                end
                end
            end
        end

        function [ip]         = get_irreducible_cluster(pc,nvertices,cutoff)
            
            import am_lib.* am_dft.*

            % construct concrete supercell for determining pairs
            [uc,uc.u2p,uc.p2u] = get_supercell(pc, diag(ceil(2*cutoff./normc_(pc.bas))) ); uc.u2i = pc.p2i(uc.u2p);
            
            % [pc-frac] get symmetries
            [~,~,S] = get_symmetries(pc); Q = get_cluster_symmetries(S,nvertices);

            % get all possible clusters with natoms
            [L{1:(nvertices-1)}]=deal([1:uc.natoms]); [Y{nvertices:-1:1}] = ndgrid(L{:},uc.p2u); 
            V = reshape(cat(nvertices+1,Y{1:nvertices}),[],nvertices).';
            
            % [cart] exclude any cluster for which at least one bond length is longer than the cutoff
            d_cart_ = @(dX) normc_(uc.bas*dX); bond_ij = nchoosek_(nvertices,2); nbonds = size(bond_ij,2); ex_ = true(1,size(V,2));            
            for i = 1:nbonds; ex_(ex_) = d_cart_( uc.tau(:,V(bond_ij(1,i),ex_)) - uc.tau(:,V(bond_ij(2,i),ex_)) )<cutoff; end
            
            % [pc-frac] create cluster tau = [X, natoms, nclusters]
            X = [uc.tau2pc*uc.tau;uc.species;uc.u2i;uc.u2p]; V=V(:,ex_);

            % get supercell cluster
            sp_ = @(pc,X,V,cutoff,Q) struct('model',[],...
                'units','frac-pc','bas',pc.bas, ...
                'symb',{pc.symb},'mass',pc.mass, ...
                'x2p',X(6,:),'x2i',X(5,:),'species',X(4,:),'tau',X(1:3,:), ...
                'nQs',size(Q{1},3),'Q',{Q},'Dj',[], ...
                'cutoff',cutoff,'nvertices',size(V,1),'nclusters',size(V,2),'cluster',V);
            sp = sp_(pc, X, V, cutoff, Q);

            % get action of symmetry operations on supercell cluster
            [~,i2p] = get_action(sp);

            % build
            V = sp.cluster(:,i2p); [Z,~,IC] = uniquec_(V(:).'); V=reshape(IC,size(V)); X=X(:,Z);

            % build irreducible clusters 
            ip_ = @(pc,X,V,cutoff,Q) struct('model',[],...
                'units','frac-pc','bas',pc.bas, ...
                'symb',{pc.symb},'mass',pc.mass, ...
                'x2p',X(6,:),'x2i',X(5,:),'species',X(4,:),'tau',X(1:3,:),  ...
                'nQs',size(Q{1},3),'Q',{Q},'Dj',[], ...
                'cutoff',cutoff,'nvertices',size(V,1),'nclusters',size(V,2),'cluster',V);
            ip = ip_(pc, X, V, cutoff, Q);
            
            % now is the time to sort based on distances
            REF = reshape(ip.bas*ip.tau(:,circshift(ip.cluster,1,1)),[3*ip.nvertices,ip.nclusters]);
            POS = reshape(ip.bas*ip.tau(:,          ip.cluster     ),[3*ip.nvertices,ip.nclusters]);
            d = normc_(REF-POS)/sqrt(ip.nvertices); fwd = rankc_([d;min(sign(POS));min(sign(POS))]);
            ip.cluster = ip.cluster(:,fwd);

            function [Q,nQs]      = get_cluster_symmetries(S,nvertices)
                % combine space symmetry with permutation of atomic positions
                M = perms([nvertices:-1:1]).'; Q{1} = repmat(S,1,1,size(M,2)); Q{2} = repelem(M,1,size(S,3)); nQs=size(Q{1},3);
            end
        end

        function [pp]         = get_primitive_cluster(ip,uc)
        
            import am_lib.* am_dft.*

            % [pc-frac] create cluster tau = [X, natoms, nclusters]
            X = [ip.tau;ip.species]; tau = reshape( X(:,ip.cluster), size(X,1), ip.nvertices, ip.nclusters);

            % [pc-frac] apply transformation tau = [X, natoms, nclusters, nsymmetries]
            tau = apply_symmetry(ip.Q, tau); 
            
            % [pc-frac] shift reference atom to primitive cell 
            tau(1:3,:,:,:) = tau(1:3,:,:,:) - floor(tau(1:3,1,:,:));

            % [pc-frac] make a list of unique positions and assign each a number (rank by species then by position)
            X = uniquec_( reshape(tau,size(tau,1),[]) ); 

            % get primitive and irreducible labels for positions
            fwd = member_(mod_(X(1:3,:)),mod_(ip.tau)); x2p = ip.x2p( fwd ); x2i = ip.x2i( fwd );
            
            % assign clusters unique labels
            [V,~,V_p2i]=uniquec_( member_(tau/10,X/10) );

            % get irreducible cluster indicies by connecting symmetrically equivalent clusters with a graph
            PM = reshape(V_p2i, [ip.nclusters, ip.nQs] ); [~,ip2pp,pp2ip] = get_connectivity( PM ); 

            % get symmetries which take irrep to orbit
            pp_nclusters = max(PM(:)); i2o = zeros(1,pp_nclusters);
            for i = 1:pp_nclusters; i2o(i) = max(findrow_(PM==i)); end

            % create structure
            pp_ = @(ip,X,V,x2p,x2i,pp2ip,ip2pp,Q) struct('units','frac-pc','bas',ip.bas,...
                'symb',{ip.symb},'mass',ip.mass,...
                'x2p',x2p,'x2i',x2i,'species',X(4,:),'tau',X(1:3,:), ...
                'nQs',size(Q{1},3),'Q',{Q}, ...
                'cutoff',ip.cutoff,'nvertices',size(V,1),'nclusters',size(V,2),'cluster', V,'i2o',[],'o2i',[], ...
                'pp2ip',pp2ip,'ip2pp',ip2pp);
            pp = pp_(ip,X,V,x2p,x2i,pp2ip,ip2pp,ip.Q);

            % get inverse elements from multiplication table 
            [~,~,I] = get_multiplication_table(pp.Q); [PM,i2p,p2i]=get_action(pp);
            
            % figure out which symmetries take the orbit to the irrep
            pp.o2i = findrow_(PM==i2p(p2i).').'; % accessc_(PM2.',o2i(:).')
            pp.i2o = I(pp.o2i).';
            
            % make sure that ip2pp matches exactly (if it doesn't check SORT above)
            for i = 1:ip.nclusters
                if ~all_(eq_( ip.tau(:,ip.cluster(:,i)), pp.tau(:,pp.cluster(:,pp.ip2pp(i)))))
                    error('ip and pp clusters mismatch');
                end
            end
            % make sure first occurance of pp matches ip
            if ~all_(eq_( pp.ip2pp, find(pp.o2i==1) ))
                error('first occuranc of equivalent pp must match ip');
            end
            
            % get matrix element dimensions (to be used below)
            for matrix_element_ = {'v_tb','v_fc'}; if isfield(ip,matrix_element_{:})
                d(ip.x2p(ip.cluster(1,:))) = cellfun(@(x)size(x,1),ip.(matrix_element_{:})); 
                d(ip.x2p(ip.cluster(2,:))) = cellfun(@(x)size(x,2),ip.(matrix_element_{:}));
            end; end
        
            % expand matrix elements as well
            for matrix_element_ = {'v_tb','v_fc','fc','tb'}
                if isfield(ip,matrix_element_{:})
                    % get dimensions
                    d = cellfun(@(x)size(x,1),ip.Dj);
                    % loop over clusters
                    for p = 1:pp.nclusters
                        c = pp.pp2ip(p);
                        % construct symmetry
                        qi = pp.i2o(p); ijk = ip.x2i(ip.cluster(:,c));
                        W = 1; for j = ijk; W = am_lib.kron_( ip.Dj{j}(:,:,qi), W ); end 
                        W = am_lib.transpose_( d(ijk), ip.Q{2}(:,qi) ) * W; W = am_lib.wdv_(W);
                        % get irreducible matrices
                        if contains(matrix_element_{:},'v') % symbolic
                            M = ip.(matrix_element_{:}){c}; W = wdv_(W);
                        else                                % numeric
                            M = double(ip.W{c}*ip.fc{c}(:));
                        end
                        % transform to primitive matrices
                        pp.(matrix_element_{:}){p} = reshape( W * M(:) , d(ijk(ip.Q{2}(:,qi))) );
                    end
                end
            end

            % if uc is present, link pp to uc
            if nargin > 1
                pp.ncenters=[];pp.norbits=[]; 
                [pp.c_id,pp.o_id,pp.q_id,pp.pp_id,pp.ip_id] = get_cluster_map(pp,uc);
                [pp.norbits,pp.ncenters] = cellfun(@(x)size(x),pp.o_id);
            end

            function [c_id,o_id,q_id,pp_id,ip_id] = get_cluster_map(pp,uc)
                % convert positions to [uc-frac]
                pp_tau = uc.bas\pp.bas*pp.tau;
                pc_tau = uc.tau(:,uc.p2u);
                uc_tau = uc.tau;
                pc2pp  = am_lib.member_(pc_tau,pp_tau); % indicies of pc atoms in pp
                pp_id  = [];

                % loop over primitive cell atoms
                for m = 1:numel(uc.p2u)
                    % record unit cell atoms of primitive type m
                    c_id{m} = find(uc.u2p==m); ncenters = numel(c_id{m});
                    % find orbits around c_id{m}(n), count their numbers
                    ex_ = [pp.cluster(1,:)==pc2pp(m)]; norbits = sum(ex_); 
                    % save ip and pp indicies
                    pp_id{m} = find(ex_(:));
                    ip_id{m} = pp.pp2ip(ex_).';
                    % allocate space
                    o_id{m} = zeros(norbits,ncenters,pp.nvertices-1);
                    q_id{m}(1:norbits,1) =  pp.o2i(ex_);
                    % loop over centers
                    for n = 1:ncenters
                        % center clusters atoms on uc reference frame
                        pp_tau = pp_tau - pp_tau(:,pc2pp(m)) + uc_tau(:,c_id{m}(n));
                        % shift to positive octant
                        pp_tau = am_lib.mod_(pp_tau);
                        % compare to get indicies of pp atoms in uc
                        o_id{m}(:,n,:) = reshape( am_lib.member_(pp_tau(:,pp.cluster(2:end,ex_).'),uc_tau), [norbits , pp.nvertices-1]);
                        % o_id{m}(:,n,1) =  am_lib.member_(pp_tau(:,pp.cluster(2,ex_)),uc_tau);
                        % o_id{m}(:,n,2) =  am_lib.member_(pp_tau(:,pp.cluster(2,ex_)),uc_tau);
                    end
                end
            end
        end

        function [PM,i2p,p2i,s_ck,g_ck] = get_action(ip)
            import am_dft.* am_lib.*
            
            % get field properly
            if     isfield(ip,'S'); sym='S'; nsyms='nSs';
            elseif isfield(ip,'Q'); sym='Q'; nsyms='nQs'; end
            % composite vector
            X = [ip.tau;ip.species]; Xd = size(X,1);
            % apply symmetry
            X_action = apply_symmetry(ip.(sym), reshape(X(:,ip.cluster), Xd, ip.nvertices, ip.nclusters) );
            % shift first atom to primitive
            X_action(1:3,:,:,:) = X_action(1:3,:,:,:) - floor(X_action(1:3,1,:,:));
            % reassign index
            [~,~,V_p2i] = uniquec_( reshape(X_action, Xd*ip.nvertices, ip.nclusters*ip.(nsyms)) );
            % get permutation matrix
            PM = reshape(V_p2i, ip.nclusters, ip.(nsyms));
            % get connectivity
            [~,i2p,p2i] = get_connectivity( PM );
            % get stabilzier and generators [nQs, nclusters]
            s_ck = false(ip.(nsyms),ip.nclusters);
            g_ck = false(ip.(nsyms),ip.nclusters);
            s_ck = [PM(i2p,:)==PM(i2p,1)].';
            for i = 1:numel(i2p); [~,a,~]=unique(PM(i2p(i),:)); g_ck(a,i)=true; end
        end

        function [bz]         = get_dispersion(ip,bz)
            % Get phonon and electron dispersions on bz.
            % model is either bvk or tb
            %
            % Q: Are the eigenvalues and eigenvectors the same at
            %    symmetry-equivalent positions in reciprocal space?
            % A: Eigenvalues, yes; eigenvectors are (to within a sign)
            %    rotated by the same point symmetry relating the different
            %    k-points. Check it out with:
            %
            %     % generate point symmetries [recp-frac]
            %     [~,~,~,R] = get_symmetries(pc); nRs = size(R,3); R = permute(R,[2,1,3]);
            %     % get orbit of a rank k point
            %     k = matmul_(R,rand(3,1)); nks = size(k,2);
            %     % get eigenvectors and eigenvalues at each orbit point
            %     for i = 1:nks
            %     input = num2cell([bvk.fc{:},(pc.bas.'\k(:,i)).',bvk.mass]); % [cart]
            %     [U(:,:,i),hw(:,i)] = eig( force_hermiticity_(bvk.D(input{:})) ,'vector');
            %     end
            %
            %     % covert symmetries [recp-frac] -> [recp-cart]
            %     sym_rebase_ = @(B,R) matmul_(matmul_(B,R),inv(B));
            %     R_cart = sym_rebase_(inv(pc.bas).',R);
            %
            %     % rotate the first eigenvector and compare with other kpoints
            %     X = matmul_( blkdiag_(R_cart,R_cart) ,  U(:,1,1) );
            %     abs(X) - abs(squeeze(U(:,1,:))) < am_lib.eps
            %
            %     % compare eigenvalues
            %     hw-hw(:,1) < am_lib.eps
            %

            import am_lib.* am_dft.*

            fprintf(' ... computing dispersion '); tic;
            if     contains(ip.model,'tb')
                % get eigenvalues
                bz.nbands = ip.nbands;
                bz.E = zeros(bz.nbands,bz.nks);
                bz.V = zeros(bz.nbands,bz.nbands,bz.nks);
                for i = 1:bz.nks
                    % define input ...
                    input = num2cell([ip.vsk{:},bz.k(:,i).']);
                    % ... and evaluate (V are column vectors)
                    [bz.V(:,:,i),bz.E(:,i)] = eig(  ip.H(input{:}) ,'vector');
                    % check hermicity of hamiltonian
                    if any(imag(bz.E(:,i)))>am_dft.eps; fprintf('\n'); error('H is not Hermitian (%g,%i)\n',max(imag(bz.E(:,i))),i); end
                    % sort energies
                    [~,fwd]=sort(bz.E(:,i)); bz.E(:,i)=bz.E(fwd,i); bz.V(:,:,i)=bz.V(:,fwd,i);
                end
            elseif contains(ip.model,'bvk')
                % get eigenvalues
                bz.nbranches = ip.nbranches;
                bz.hw = zeros(bz.nbranches,bz.nks);
                bz.U  = zeros(bz.nbranches,bz.nbranches,bz.nks);
                for i = 1:bz.nks
                    % define input ...
                    % input = num2cell([bvk.fc{:},bz.k(:,i).',bvk.mass]); % [pc-frac]
                    input = num2cell([ip.fc{:},bz.k(:,i).',ip.mass]); % [cart]
                    % ... and evaluate (U are column vectors)
                    [bz.U(:,:,i),bz.hw(:,i)] = eig( ip.D(input{:}) ,'vector');
                    % check hermicity of hamiltonian
                    if any(imag(bz.hw(:,i))>am_dft.eps); fprintf('\n'); error('H is not Hermitian (%e,%i)\n',max(imag(bz.hw(:,i))),i); end
                    % sort energies
                    [~,fwd]=sort(bz.hw(:,i)); bz.hw(:,i)=bz.hw(fwd,i); bz.U(:,:,i)=bz.U(:,fwd,i);
                    % correct units
                    bz.hw(:,i) = sqrt(real(bz.hw(:,i))) * am_lib.units_eV;
                end
            end
            fprintf('(%.f secs)\n',toc);
        end

        function                plot_dispersion(ip, bzp, varargin)
            % model is either bvk or tb

            import am_lib.* am_dft.*

            % define figure properties
            fig_ = @(h)       set(h,'color','white');
            axs_ = @(h,qt,ql) set(h,'Box','on','XTick',qt,'Xticklabel',ql);
            fig_(gcf);

            if     contains(ip.model,'tb')
                % get electron band structure along path
                bzp = get_dispersion(ip,bzp);
                % plot results
                plot(bzp.x,sort(bzp.E), varargin{:});
                axs_(gca,bzp.qt,bzp.ql); axis tight; ylabel('Energy [eV]'); xlabel('Wavevector k');
            elseif contains(ip.model,'bvk')
                % get phonon band structure along path
                bzp = get_dispersion(ip,bzp);
                % plot results
                plot(bzp.x,sort(real(bzp.hw)*1E3),'-k',bzp.x,-sort(abs(imag(bzp.hw))), varargin{:});
                axs_(gca,bzp.qt,bzp.ql); axis tight; ylabel('Energy [meV]'); xlabel('Wavevector k');
            else 
                error('model unknown');
            end
        end


        % phonons (harmonic)

        function [md]         = run_bvk_md(bvk,pp,uc,dt,nsteps,Q,T)
            % set time step [ps ~ 0.1], number of MDs steps, Nose-Hoover "mass" Q, and temperature T [K]
            % dt = 0.1; nsteps = 10000; Q = 1; T = 300; [md] = run_bvk_md(bvk,pp,uc,dt,nsteps,Q,T)

            import am_lib.*
            import am_dft.*

            % build force constants
            for m = 1:pp.pc_natoms
                phi{m} = zeros(3,3*pp.npairs(m));
            for j = 1:pp.npairs(m)
                % get indicies
                i = pp.i{m}(j); iq = pp.iq{m}(j);
                % get irrep force constant indicies
                iphi = reshape(bvk.W{i}*bvk.fc{i}(:),3,3);
                % rotate force constants from irrep to orbit
                phi{m}(1:3,[1:3]+3*(j-1)) = pp.Q{1}(1:3,1:3,iq) * permute(iphi,pp.Q{2}(:,iq)) * pp.Q{1}(1:3,1:3,iq).';
            end
            end

            % allocate arrays and constants [cart]
            k_boltz = 8.6173303E-5; % [eV/K]
            u = single(zeros(3,uc.natoms,nsteps)); KE=zeros(1,nsteps);
            v = single(zeros(3,uc.natoms,nsteps)); PE=zeros(1,nsteps);
            f = single(zeros(3,uc.natoms,nsteps));
            a = single(zeros(3,uc.natoms,nsteps));

            % set initial value conditions: small displacement just to get atoms moving
            u(:,:,1) = (.5-rand(3,uc.natoms))*0.00001;
            v(:,:,1) =    zeros(3,uc.natoms);

            % % set velocity gaussianly distributed sqrt( k_B * K / amu) = 0.000911856 [Ang/fs]
            % r1 = rand(3,uc.natoms); r2 = rand(3,uc.natoms); u(:,:,1) = zeros(3,uc.natoms);
            % v(:,:,1) = sqrt(-2*log(r1)).*cos(2*pi*r2).*sqrt(T./uc.mass(uc.species)/2) * 0.000911856;

            % run md using verlet algorithm
            fprintf('%10s   %10s     %10s   %10s   %10s \n','step','temp','PE','KE','PE+KE');
            for j = 1:nsteps
                % compute force [eV/Ang]
                for m = 1:pp.pc_natoms; f(:,pp.c{m},j) = - phi{m} * reshape(u(:,pp.o{m},j), size(pp.o{m}).*[3,1]); end

                % compute potential energy [eV/Ang * Ang -> eV]
                PE(j) = - flatten_(u(:,:,j)).'*flatten_(f(:,:,j));

                % 4) compute kinetic energy : amu * (Ang/fs)^2 -> 103.6382 eV
                KE(j) = uc.mass(uc.species)*sum((uc.bas*v(:,:,j)).^2,1).'/2 * 103.6382;

                % 6) current temperature
                Tj = 2/3*KE(j)/uc.natoms/k_boltz;

                % 5) compute Nose-Hoover drag: p_eta = KE - TE
                nosehoover = v(:,:,j)/Q * ( Tj - T ) / uc.natoms;

                % 6) get acceleration
                a(:,:,j+1) = f(:,:,j) ./ uc.mass(uc.species);

                % ***) update md [frac]: x' = x + v * dt; v' = v + a * dt; Nose-Hoover dv/dt becomes a - p_eta / Q * v;
                if j ~= nsteps
                    u(:,:,j+1) = u(:,:,j) + v(:,:,j) .* dt + a(:,:,j)./2 .* dt^2;
                    v(:,:,j+1) = v(:,:,j) + dt * ( (a(:,:,j)+a(:,:,j+1))./2 - nosehoover);
                end

                % print
                if mod(j,50)==1
                    fprintf('%10i   %10f K   %10f   %10f   %10f \n',j,Tj,PE(j),KE(j),PE(j)+KE(j));
                    %
                    figure(1); set(gcf,'color','white');
                    plot([1:nsteps],[KE;PE;KE+PE].');legend('KE','PE','KE+PE')
                    xlabel('step'); ylabel('energy [eV]'); xlim([0,nsteps]); drawnow;
                end
            end

            % convert [cart] to [frac] and u to tau
            tau = matmul_(inv(uc.bas),u)+uc.tau;
            v   = matmul_(inv(uc.bas),v);
            f   = matmul_(inv(uc.bas),f);

            % define md creation function [frac]
            md_ = @(uc,force,tau,vel,dt) struct('units','tau=frac-recp; bas=ang',...
                'bas',uc.bas,'symb',{{uc.symb{:}}},'mass',uc.mass,'nspecies',uc.nspecies, ...
                'natoms',uc.natoms,'force',force,'tau',tau,'vel',vel,'species',uc.species, ...
                'dt',dt,'nsteps',size(tau,3));
            md = md_(uc,f,tau,v,dt);
        end

        function [bvk]        = interpolate_bvk(bvk_1,bvk_2,n)
            % interpolates force constants and masses from bvk_1 and bvk_2 on n points (includes end points)

            import am_lib.* am_dft.*

            bvk_ = @(bvk,mass,fc) struct('units','tau=cart; bas=ang','bas',bvk.bas,'recbas',bvk.recbas,'natoms',bvk.natoms,'mass',mass, ...
                'nshells',bvk.nshells,'W',{bvk.W},'shell',{bvk.shell},'nbranches',bvk.nbranches,'D',bvk.D,'fc',{fc});

            fc_interp = nlinspace( [bvk_1.fc{:}] , [bvk_2.fc{:}] , n );
            mu_interp = nlinspace( [bvk_1.mass]  , [bvk_2.mass]  , n );

            % get dimensions
            for i = 1:numel(bvk_1.fc); m(i) = numel(bvk_1.fc{i}); end; E = cumsum(m); S = E - m + 1;

            % create bvks cells
            for i = 1:n
                for j = 1:numel(E); fc{j} = fc_interp(S(j):E(j),i).'; end
                bvk(i) = bvk_(bvk_1,mu_interp(:,i).',fc);
            end
        end

        function                plot_bvk_vs_aimd(uc,md,ip,bvt,pt)

            import am_lib.* am_dft.*

            % select algo 2 which is faster
            algo = 2;
            
            % get primitive pairs
            pp = get_primitive_cluster(ip,uc);
            % [cart] get displacements
            u = matmul_( md.bas, mod_( md.tau-uc.tau +.5 )-.5 );
            % [cart] get aimd forces
            f_aimd= matmul_( md.bas, md.force );
            % [cart]  get harmonic forces
            f_har = get_bvk_forces(ip, uc, u, algo);
            % [cart] get anharmnoic forces
            f_anh = 0;
            if nargin == 6; if and(~isempty(bvt),~isempty(pt))
            f_anh = get_bvt_forces(bvt,pt,u,algo);
            end; end

            % plot correlation
            plotcorr_(f_anh(:)+f_har(:),f_aimd(:)); xlabel('aimd force [eV/Ang]'); ylabel('bvk force [eV/Ang]');
        end

        function [T,KE,PE,q_sk] = plot_md_stats(uc,md,fbz,bvk)

            import am_lib.* am_dft.*

            % match uc to md
            uc = match_cell(uc,md);

            % [cart] get displacements and forces; [amu * (Ang/fs)^2 ] = 103.6382 [eV]
            u = matmul_( md.bas, mod_( md.tau-uc.tau +.5 )-.5 );
            f = matmul_( md.bas, md.force );
            v = matmul_( md.bas, md.vel ) * sqrt(103.6382);

            % get potentoal energy : [eV/Ang * Ang] = [eV]
            PE(:,1) = -dot(reshape(u,[],md.nsteps),reshape(f,[],md.nsteps),1)/2;
            % get kinetic energy : [amu * (Ang/fs)^2 ] = 103.6382 [eV] (incorporate above in to v)
            KE(:,1) =  reshape(sum(uc.mass(uc.species).*dot(v,v,1),2),1,[])/2;

            % get temperature : amu * (Ang/fs)^2/ k_B = 1.20267E6 K
            k_boltz = 8.6173303E-5; T = 2/3*KE/uc.natoms/k_boltz;

            % perform normal-mode analysis
            if nargin == 4
                % convert phonon energies back to wierd units
                hw = real(fbz.hw)./am_lib.units_eV;

                % get normal transformations from uc eigenvector (Wallace p 115 eq 10.48)
                U   = expand_bvk_eigenvectors(bvk,uc,fbz);
                u2q = U'  .* repelem(sqrt(uc.mass(uc.species)).',3,1).';
                v2p = U.' .* repelem(sqrt(uc.mass(uc.species)).',3,1).';

                % get normal modes (Wallace p 115 eq 10.49)
                q_sk = (u2q*reshape(u,[],md.nsteps));
                p_sk = (v2p*reshape(v,[],md.nsteps));

                % get energies (Wallace p 115 eq 10.53)
                PE(:,2) = dot(abs(q_sk).*hw(:),abs(q_sk).*hw(:),1)/2;
                KE(:,2) = dot(abs(p_sk)       ,abs(p_sk)       ,1)/2;
            end

            % if no output is requested, plot and print stuff
            if nargout == 0
                % print a few time steps
                Z = [[1:md.nsteps].',T(:,1),PE(:,1),KE(:,1),KE(:,1)+PE(:,1)];
                fprintf('%10s   %10s   %10s   %10s   %10s \n','step [#]','T [K]','PE [eV]','KE [eV]','PE+KE [eV]');
                if md.nsteps<100; fprintf('%10i   %10f   %10f   %10f   %10f \n',Z(1:end,:).'); else
                                  fprintf('%10i   %10f   %10f   %10f   %10f \n',Z(1:99,:).');
                                  fprintf('%10i   %10f   %10f   %10f   %10f \n',Z(100:50:end,:).');
                end

                % plot energies
                set(gcf,'color','w');
                if nargin == 4
                subplot(3,4,1:3); plot(1:md.nsteps,KE(:,1),'-',1:md.nsteps,PE(:,1),'-',1:md.nsteps,KE(:,1)+PE(:,1),'-',...
                                       1:md.nsteps,KE(:,2),'.',1:md.nsteps,PE(:,2),'.',1:md.nsteps,KE(:,2)+PE(:,2),'.');
                                legend('KE','PE','KE+PE','nKE','nPE','nKE+nPE');
                                xlabel('time step'); ylabel('energy [eV]'); axis tight;
                subplot(3,4,4); plot3(repelem([1:md.nsteps].',1,fbz.nks*bvk.nbranches),real(q_sk.*hw(:)).',real(p_sk).');
                                view([1 0 0]); box on; ax=axis; maxax=max(abs(ax(3:6))); axis([ax(1:2),[-1 1 -1 1].*maxax]); daspect([md.nsteps/maxax 1 1]);
                else
                subplot(3,1,1); plot(1:md.nsteps,KE(:,1),'-',1:md.nsteps,KE(:,2),'o',1:md.nsteps,KE(:,1)+PE(:,1),'^');
                                legend('KE','PE','KE+PE');
                                xlabel('time step'); ylabel('energy [eV]'); axis tight;
                end

                % plot position and velocity histograms
                nbins  = 101;
                dist_v = reshape(normc_(v),uc.natoms,md.nsteps);
                dist_u = reshape(normc_(u),uc.natoms,md.nsteps);
                bin_v  = linspace(0,max(dist_v(:)),nbins);
                bin_u  = linspace(0,max(dist_u(:)),nbins);
                hist_v = zeros(nbins-1,md.nsteps);
                hist_u = zeros(nbins-1,md.nsteps);
                for i = 1:md.nsteps
                    [hist_v(:,i)] = histcounts(dist_v(:,i),bin_v);
                    [hist_u(:,i)] = histcounts(dist_u(:,i),bin_u);
                end
                [Y{1:2}]=meshgrid(1:md.nsteps,bin_v(1:(nbins-1))); hist_v(:,1)=[]; Y{1}(:,1)=[]; Y{2}(:,1)=[];
                subplot(3,1,2); surf(Y{1},Y{2},hist_v,'Facecolor','interp','edgecolor','none');
                                view(2); axis tight; xlabel('time step'); ylabel('velocity [Ang/fs]');
                [Y{1:2}]=meshgrid(1:md.nsteps,bin_u(1:(nbins-1))); hist_u(:,1)=[]; Y{1}(:,1)=[]; Y{2}(:,1)=[];
                subplot(3,1,3); surf(Y{1},Y{2},hist_u,'Facecolor','interp','edgecolor','none');
                                view(2); axis tight; xlabel('time step'); ylabel('displacement [Ang]');

            end
        end


        % phonons (anharmonic)

        function [bvt,bvk]    = get_bvt_force_constants(uc,md,bvk,pp,bvt,pt)
            % Extracts symmetry adapted third-order force constant and
            % second order force constants
            %
            % Q: How to ensure that the extraction process is working
            %    correctly? Is it working correctly?
            % A: Yes, it works. Checked it by doing the reverse operation and
            %    extracting force constants from a known displacement-force
            %    set using:
            %
            %         u = matmul_( md.bas, mod_( md.tau-uc.tau +.5 )-.5 );
            %         md.force = matmul_( inv(md.bas), get_bvt_forces(bvt,pt,u,1) + get_bvk_forces(bvk,pp,u,1) );
            %         [cmp] = get_bvt_force_constants(uc,md,bvk,pp,bvt,pt);
            %         plot([bvt.fc{:}],[bvt2.fc{:}],'.')
            %
            % Q: How to factor out displacements from third-rank tensors?
            % A: Need to create the matrix
            %
            %       [ u_1*w_1, u_2*w_1, u_3*w_1, u_1*w_2, u_2*w_2, u_3*w_2, u_1*w_3, u_2*w_3, u_3*w_3,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0]
            %       [       0,       0,       0,       0,       0,       0,       0,       0,       0, u_1*w_1, u_2*w_1, u_3*w_1, u_1*w_2, u_2*w_2, u_3*w_2, u_1*w_3, u_2*w_3, u_3*w_3,       0,       0,       0,       0,       0,       0,       0,       0,       0]
            %       [       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0, u_1*w_1, u_2*w_1, u_3*w_1, u_1*w_2, u_2*w_2, u_3*w_2, u_1*w_3, u_2*w_3, u_3*w_3]
            %
            %     using the Z matrix below:
            %
            %       Z = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            %            0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            %            0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            %            0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            %            0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;
            %            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;
            %            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;
            %            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0;
            %            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1].';
            %
            %    That is: u = sym('u_%d',[3,1]); w = sym('w_%d',[3,1]); reshape(Z* reshape(u*w.',[],1),3,27)
            %

            import am_lib.* am_dft.*

            % [cart] get displacements and forces
            u = matmul_( md.bas, mod_( md.tau-uc.tau +.5 )-.5 );
            f = matmul_( md.bas, md.force );

            % choose the algorthim to extract third order anharmonic force
            % constants; algorithm 2 is more accurate
            switch 2
                case 1
                    % subtract harmonic forces first than get anharmonic
                    % force constants from the residuals

                    % [cart] subtract harmonic forces
                    if and(~isempty(bvk),~isempty(pp))
                        f = f - get_bvk_forces(bvk,pp,u);
                    end

                    % get U matrix and indexing I
                    [U,I] = get_bvt_U_matrix(bvt,pt,u);

                    % [cart] solve for force constants
                    fc = - U \ f(I);

                    % save force constants
                    s_id = repelem([1:bvt.nshells],cellfun(@(x)size(x,2),bvt.W));
                    for s = 1:bvt.nshells; bvt.fc{s} = double(fc(s_id==s)).'; end

                case 2
                    % get the second order harmonic and third order anharmonic
                    % force constants all in one go without first subtrating the
                    % harmonic forces.

                    % get U matrix and indexing I
                    [Uh,Ih] = get_U_matrix(bvk,pp,u);
                    % get U matrix and indexing I
                    [Ua,Ia] = get_bvt_U_matrix(bvt,pt,u);
                    % check that indicies match
                    if maxabs_(Ih-Ia)>am_lib.eps
                        error('indicies do not match');
                        % this case has not yet been coded...
                    else
                        I = Ih;
                    end

                    % solve for force constants
                    fc = - [ Uh, Ua ] \ f(I);

                    % save harmonic and anharmnoic force constants
                    h_id = repelem([1:bvk.nshells],cellfun(@(x)size(x,2),bvk.W));
                    a_id = repelem([1:bvt.nshells],cellfun(@(x)size(x,2),bvt.W));
                    s_id = [h_id,bvk.nshells+a_id];
                    for s = 1:bvk.nshells; bvk.fc{s} = double(fc(s_id==s            ).'); end
                    for s = 1:bvt.nshells; bvt.fc{s} = double(fc(s_id==s+bvk.nshells).'); end
            end
        end


        % tight binding
  
        function                plot_tb_vs_dft(tb,dft)

            import am_lib.* am_dft.*

            % plot correlation
            plotcorr_( flatten_( dft.E([1:tb.nbands]+tb.nskips,:)        ) , ...
                       flatten_( eval_energies_(tb,[tb.vsk{:}],dft.k) ) );
            xlabel('dft energies [eV]'); ylabel('tb energies [eV]');
        end

        
        % tensors
        
        function [M] = get_second_rank_tensor_symmetry(R)
            % get sizes
            n = size(R,1); m = size(R,1); d = [n,m];
            % build equations
            W = am_lib.kron_( R , R ); 
            W = sum( W - eye(n^2,m^2), 3);
            % get nullspace
            W = null(W); W = am_lib.frref_(W.').'; W = W./am_lib.accessc_(W,am_lib.findrow_(W.').'); W = am_lib.wdv_(W); 
            % % define parameters
            C = sym('c%d%d', d); C = C(am_lib.findrow_(W.'));
            % get tensor 
            M = reshape( sym(W)*C(:), d );
        end
    end

    
    % diffraction
    
    methods (Static)
        
        function [I]       = simulate_xrr_engine_(eta,th,hv,thickness,filling,roughness,sample_length,algo)
            %
            % R = simulate_xray_reflectivity(layer,th,hv,thickness,filling,roughness)
            % 
            % eta       [unit-less]     x-ray refractive index
            % th        [deg]           angles
            % thickness [nm]            thickness
            % filling   [unitless]      multiplies density
            % roughness [nm]            inteface roughness
            % method    [1,2]           transfer matrix (explicit,slow), recursive parratt (default,fast)
            
            import am_lib.*
            
            if nargin < 8; algo = 2; end
            if nargin < 7; sample_length = []; end
            if size(th,2)>size(th,1); th=th.'; end 
            
            % get number of layers and number of angles
            nlayers = numel(thickness); nths = numel(th); kz = zeros(nths,nlayers);
            
            % get photon wavelength
            lambda  = get_photon_wavelength(hv);
            
            % get kx and kz
            % Note: k are incident wavevectors and not diffraction vectors
            get_kz = @(th,lambda) sind(th)./lambda;
            get_kx = @(th,lambda) cosd(th)./lambda;

            switch algo
                case 1 % explicit transfer matrix (slower)
                    % [nths,nlayers] solve wavevector boundary conditions to get
                    % out-of-plane wavevector component in layer kz [nths,nlayers]
                    % 1. k2  = (n2/n1) k1  is  snell's law
                    % 2. k1x = k2x
                    get_kz_in_layer = @(n1,n2,k1z,k1x) sqrt( (n2./n1).^2 .* k1z.^2 + ((n2./n1).^2-1) .* k1x.^2 );
                    for i = 1:nlayers
                        kz(:,i) = get_kz_in_layer(1,filling(i).*(eta(:,i)-1)+1,get_kz(th,lambda),get_kx(th,lambda));
                    end

                    % get transfer matricies R and T which describe propagation
                    % across interfaces and media
                    % syms a b c d k1 k2 z
                    % % boundary conditions for plane waves at an interface
                    % eq(1) = a * exp(1i*k1*z) + b * exp(-1i*k1*z) == c*exp(1i*k2*z) + d*exp(-1i*k2*z);
                    % eq(2) = diff(eq(1),z);
                    % solution = solve(eq,a,b);
                    % % transfer matrix
                    % T(1,1:2) = simplify(equationsToMatrix(solution.a,c,d));
                    % T(2,1:2) = simplify(equationsToMatrix(solution.b,c,d));
                    % % set interface at 0 for simplicity
                    % T = subs(T,z,0)
                    % simplify
                    % T = simplify(expand(T));
                    % % This is exact when there is no roughness.
                    % get_transfer_matrix_at_interface = @(k1,k2,sigma) [ ... 
                    %       [ exp(1i*( - k1 + k2 )) * (k1 + k2), exp(1i*( - k1 - k2 ))*(k1 - k2)]
                    %       [ exp(1i*( + k1 + k2 )) * (k1 - k2), exp(1i*( + k1 - k2 ))*(k1 + k2)] ] ./ (2*k1)];
                    % to add roughness, see below.
                    % hmm .. these two methods seem to be equivalent (when sigma == 1 and z == 1)? as they should be. 

                    get_transfer_matrix_at_interface = @(k1,k2,sigma) [ ...
                        [ (exp(-(k1 - k2).^2*(sigma*2*pi).^2/2) * (k1 + k2))/(2*k1), (exp(-(k1 + k2).^2*(sigma*2*pi).^2/2) * (k1 - k2))/(2*k1)]
                        [ (exp(-(k1 + k2).^2*(sigma*2*pi).^2/2) * (k1 - k2))/(2*k1), (exp(-(k1 - k2).^2*(sigma*2*pi).^2/2) * (k1 + k2))/(2*k1)] ];

                    get_transfer_matrix_in_medium = @(k,l) [ ...
                            [ exp(-2i.*pi.*k.*l), 0]
                            [ 0, exp(+2i.*pi.*k.*l)]];

                    % get reflection at each theta value
                    I = zeros(nths,1);
                    for j = 1:nths
                        % vacuum/first layer
                        M{1} =        get_transfer_matrix_at_interface(get_kz(th(j),lambda), kz(j,1), roughness(1) );
                        M{1} = M{1} * get_transfer_matrix_in_medium(                         kz(j,1), thickness(1) );
                        % first layer ... nth layer
                        if nlayers > 1
                            for i = 2:nlayers
                                M{i} =        get_transfer_matrix_at_interface( kz(j,i-1)  , kz(j,i), roughness(i) );
                                M{i} = M{i} * get_transfer_matrix_in_medium(                 kz(j,i), thickness(i) );
                            end
                        end
                        % apply 
                        TM = mtimes_(M{:});
                        % calculate intensity
                        I(j) = abs(TM(2,1)./TM(1,1)).^2;
                    end
                case 2 % recursive Parratt method without transfer matrix (faster)
                    thickness = thickness*2;
                    roughness = roughness*4*pi;
                    % initialize matrices
                    r = zeros(nths,nlayers); R = zeros(nths,nlayers-1);
                    %----- Wavevector transfer in each layer
                    kz(:,1) = get_kz(th,lambda);
                    for j=1:nlayers
                        kz(:,j+1)= sqrt( kz(:,1).^2 + 2*(1/lambda).^2 * (eta(:,j)-1) * filling(j) );
                    end
                    %----- Reflection coefficients (no multiple scattering)
                    for j=1:nlayers
                        r(:,j)=(  (kz(:,j)-kz(:,j+1))./(kz(:,j)+kz(:,j+1)) ) .* exp(-0.5*(kz(:,j).*kz(:,j+1))*roughness(j)^2);
                    end
                    %----- Reflectivity
                    if nlayers>1
                        R(:,1) =  (r(:,nlayers-1)  + r(:,nlayers) .* exp(2i*pi*kz(:,nlayers)*thickness(nlayers-1)) ) ...
                              ./(1+r(:,nlayers-1) .* r(:,nlayers) .* exp(2i*pi*kz(:,nlayers)*thickness(nlayers-1)) );
                    end
                    if nlayers>2; for j=2:nlayers-1
                        R(:,j) =  (r(:,nlayers-j)  + R(:,j-1) .* exp(2i*pi*kz(:,nlayers-j+1)*thickness(nlayers-j)) ) ...
                              ./(1+r(:,nlayers-j) .* R(:,j-1) .* exp(2i*pi*kz(:,nlayers-j+1)*thickness(nlayers-j)) );
                    end; end
                    %------ Intensity reflectivity
                    if nlayers==1; I = abs(r(:,1)).^2; else; I = abs(R(:,nlayers-1)).^2; end
                case 3 % vectorized
                    % Note that there are many typographical mistakes in Matt's paper.
                    thickness=thickness*2*pi; roughness=roughness*2*pi;
                    % get wave vector transfer (simplifies to: sind(th)/lambda for eta = 1); Wormington Eq 4.3
                    kz_ = @(lambda,eta,th) sqrt(eta.^2 - cosd(th).^2)/lambda;
                    % Frensel reflection coefficients; Wormington Eq 4.4
                    r_ = @(kz,roughness) (kz(:,1:(end-1))-kz(:,2:end))./(kz(:,1:(end-1))+kz(:,2:end)) .* exp(-2.*kz(:,1:(end-1)).*kz(:,2:end).*roughness.^2);
                    % Parratt recursion; Wormington Eq 4.2
                    X_ = @(X,r,kz,thickness) (r + X .* exp(2i*kz*thickness))./(1 + r.*X .* exp(2i*kz*thickness));

                    % compute wave vector and Frensel coefficients
                    kz = kz_(lambda,[ones(nths,1),(eta-1).*filling+1],th); r = r_(kz,roughness); 

                    X = r(:,nlayers);
                    for j = nlayers-1:-1:1
                        X = X_(X,r(:,j),kz(:,j+1),thickness(j));
                    end
                    
                    I = abs(X).^2;
            end
            
            % add intensity correction due to finite sample width
            if ~isempty(sample_length) && ~eq_(sample_length,0)
                xray_beam_height = 0.1; % [mm] height of x-ray beam
                th_b = asind(xray_beam_height/sample_length);
                I(th<th_b) = sind(th(th<th_b)) ./ sind(th_b) .* I(th<th_b);
            end
        end

        function [I]       = simulate_xrr(layer,th,hv,W,Imax,I0)
            import am_mbe.* am_lib.* am_dft.*
            
            if nargin<2; th = [0:0.005:4].'; end
            if nargin<3; hv = get_atomic_emission_line_energy(am_dft.get_atomic_number('Cu'),'kalpha1'); end
            if nargin<4; W=1; Imax=1; I0=0; end
            
            % get number of layers
            nlayers = numel(layer);

            % get dielectric contributions
            for i = 1:nlayers; xray_refractive_index(:,i) = get_xray_refractive_index(layer(i), hv, get_kz(th,hv)); end
            
            % compute intensity
            I = simulate_xrr_engine_(...
                xray_refractive_index, th, hv,...
                am_lib.field2array_(layer,'thickness'), ...
                am_lib.field2array_(layer,'filling')  , ...
                am_lib.field2array_(layer,'roughness'), ...
                W)*Imax+I0;
            
            % plot final result            
            semilogy(th,I,'k.-'); set(gca,'linewidth',1);set(gcf,'color','w');
            axis tight; xlabel('\theta [deg.]'); ylabel('Intensity [a.u.]');
        end

        function [I]       = simulate_xrd(layer,bz,hv,thickness,strain,flag)
            % dynamical simulation of x-ray diffraction intensity
            % i = 0; t = [30, 1E8];
            % i=i+1;layer(i) = load_cell('material','GaAs'); 
            % i=i+1;layer(i) = load_cell('material','Si');
            % bz = get_angles(layer(1),hv,100000,30,180);
            % I  = simulate_xrd(layer,bz,hv,t);
            
            import am_dft.* am_lib.* 
            
            if nargin<6; flag=''; end
            
            % get lambda an number of layers
            lambda = get_photon_wavelength(hv); nlayers = numel(layer); 
            
            % define equations according to:
            % M. Wormington, C. Panaccione, K. M. Matney, and D. K. Bowen, Philosophical Transactions of 
            % the Royal Society A: Mathematical, Physical and Engineering Sciences 357, 2827 (1999).
                % Eq 5.3, Authier Eq 8
                get_susceptibility_ = @(lambda,vol,F) (-am_lib.r_0*lambda.^2./(pi.*vol)).*F;
                % Eq 5.1 asymmetry factor
                get_asymmetry_factor_ = @(w,th2) -sind(w)./sind(th2-w);
                % Eq 5.1 polarization factor, Bartels after Eq 5
                get_polarization_ = @(th2_bragg) abs(cosd(th2_bragg));
                % Eq 5.2 angular deviation parameter, Bartels Eq 3
                get_angular_deviation_ = @(th2,th2_bragg) -4*(sind(th2/2)-sind(th2_bragg/2)).*sind(th2_bragg/2);
                % Eq 5.1 get complex deviation parameter eta, Authier Eq 25
                get_complex_deviation_ = @(alpha,b,C,chi_0,chi_H) (    alpha - chi_0.*(1-b)    ) ./ ( 2.*abs(C).*sqrt(abs(b)).*abs(chi_H) );
                % Eq 5.7; Bartels Eq 4 [OK]
                get_T_ = @(lambda, C, w, th2, chi_H, t) pi*C.*abs(chi_H)./sqrt(abs(sind(w).*sind(th2-w))) .* (t/lambda);
                % Eq 5.6; Bartels Eqs 7,8; Birkholz Eq 7.45 [OK]
                get_S1_= @(X,T,eta) (X - eta + sqrt(eta.^2-1)).*exp(-1i.*T.*sqrt(eta.^2-1));
                get_S2_= @(X,T,eta) (X - eta - sqrt(eta.^2-1)).*exp(+1i.*T.*sqrt(eta.^2-1));
                % Eq 5.5; Bartels Eq 6,9; Birkholz Eq 7.44, 7.46 Darwin-Prins formula [OK]
                get_X_infinite = @(eta)       eta -  sign(real(eta)) .* sqrt(eta.^2-1);
                get_X_finite   = @(eta,S1,S2) eta + (S1+S2)./(S1-S2) .* sqrt(eta.^2-1);
                % Bartels Eq 12, 13
                kinematical_algorithm = 1;
                switch kinematical_algorithm
                case 1
                    get_X_kinematical = @(X,eta,T)  X.*exp(-2i.*eta.*T) + (1-exp(-2i.*eta.*T))./(2.*eta);
                case 2
                    get_X_kinematical = @(X,eta,T)  X.*exp(-2i.*eta.*T) + exp(-1i.*eta.*T).*sin(eta.*T)./eta;
                end

            % evaluate recusively
            for i = nlayers:-1:1
                % get bragg points
                k = linspacen_([0;0;1],[0;0;10],10);
                k_cart_magnitude = normc_(inv(layer(i).bas .* [1 0 0; 0 1 0; 0 0 1+strain(i)] ).'*k);
                [~,th2_bragg]=kxkz2angle(0,k_cart_magnitude,hv,'w2th');
                th2_bragg = reshape(th2_bragg(th2_bragg<180),1,1,[]); 

                % consider all hkl reflections or select one
                switch 'all'
                case 'custom'
                    % SnO(001)/Al2O3(1 -1 0 2)
                    switch i
                    case 1
                        F_H   = get_structure_factor(layer(i), [0 0 1].' , hv); % SnO (0 0 1)
                    case 2
                        F_H   = get_structure_factor(layer(i), [1 -1 2].', hv); % Al2O3 (1 -1 0 2)
                    end
                case 'all'
                    F_H   = get_structure_factor(layer(i), bz    , hv); 
                case 'bede' % used in bede
                    reflection_ = 4; k = k(:,reflection_); th2_bragg=th2_bragg(reflection_);
                    F_H   = get_structure_factor(layer(i), k     , hv); 
                end

                % evaluate everything
                vol   = det(layer(i).bas .* [1 0 0; 0 1 0; 0 0 1+strain(i)]);
                C     = get_polarization_(0);
                F_0   = get_structure_factor(layer(i),[0;0;0], hv);
                chi_0 = get_susceptibility_(lambda, vol, F_0);
                chi_H = get_susceptibility_(lambda, vol, F_H);
                b     = get_asymmetry_factor_(bz.x/2, bz.x);
                alpha = get_angular_deviation_(bz.x, th2_bragg);
                eta   = get_complex_deviation_(alpha, b, C, chi_0, chi_H);

                if thickness(i) > 1E5
                    X = get_X_infinite(eta);
                else
                    if i==nlayers
                        X = zeros(1,bz.nks);
                    end
                    T  = get_T_(lambda, C, th2_bragg/2, th2_bragg, chi_H, thickness(i));
                    
                    % ad-hoc correction to make sure sizes are the same
                    % this is not the correct way to do it. ideally need to figure out how ot match each theta with each other.
                    % for this reason Taupin-Tagaki equations are not that great...
                        m = min([size(X,3),size(T,3),size(eta,3)]); 
                        X=X(:,:,1:m); T=T(:,:,1:m); eta=eta(:,:,1:m); 
                        
                    S1 = get_S1_(X,T,eta); S2 = get_S2_(X,T,eta);
                    % "kinematical" simulation
                    
                    if contains(flag,'kinematical')
                        X  = get_X_kinematical(X,eta,T);
                    else
                        X  = get_X_finite(eta,S1,S2);
                    end
                end
            end
            I = sum(abs(X).^2,3);
        end

    end

    % aux library

    methods (Static)

        % aux phonons (harmonic)
        
        function [f]     = get_bvk_forces(ip, uc, u, algo)
            %
            % Get forces f from the displacement u.
            %
            % Selecting a METHOD is optional (defaults to METHOD 2):
            %
            % METHOD 1: Build U matrix and compute forces using:
            %           F [3m * 1] = - U [3m * nFcs] * FC [ nFCs * 1] for n pairs, m atoms
            %
            % METHOD 2: Compute F directly from tensor double dot products using:
            %           F [3 * m] = - FC [3 * 3n] * U [3n * m]:
            %           for n pairs, m atoms (3 spans displacements)
            %
            % Get [cart] displacement using, for example:
            %     u = matmul_( md.bas, mod_( md.tau-uc.tau +.5 )-.5 );
            %
            % Compare forces computed using the two methos:
            % [pp]         = get_primitive_cluster(ip,uc)
            % u=rand(3,250,20);
            % tic; f1 = am_dft.get_bvk_forces(ip,pp,u,1); toc
            % tic; f2 = am_dft.get_bvk_forces(ip,pp,u,2); toc
            % plot(f1(:),f2(:),'.')
            %

            import am_lib.* am_dft.*

            if nargin ~= 4; algo=1; end
            
            switch algo
                case 1
                    % get U matrix and indexing
                    [U,I] = get_U_matrix(ip, uc, u);
                    % solve for forces f [3m * 1] = - U [3m * nFcs] * FC [ nFCs * 1]: n pairs, m atoms
                    f(I) = - U * [ip.fc{:}].';
                    % rearrange forces
                    f = reshape(f,size(u));
                case 2
                    % get sizes
                    [~,natoms,nsteps] = size(u);
                    % expand force constants over primitive pairs and link primitive pairs to unit cell
                    pp = get_primitive_cluster(ip,uc);
                    % build force constants
                    for m = 1:max(pp.x2p(:)); phi{m} = cat(2,pp.fc{ pp.pp_id{m} }); end
                    % compute forces on every atom at every step
                    f = zeros(3,natoms,nsteps);
                    for j = 1:nsteps; for m = 1:numel(pp.c_id)
                        f(1:3,pp.c_id{m},j) = - phi{m} * reshape(u(:,pp.o_id{m},j), size(pp.o_id{m}).*[3,1]);
                    end; end
            end
        end

        function [U,I]   = get_U_matrix(ip, uc, u)
            %     Formalism
            %
            %     % get Z matrix
            %     FC=sym('FC',[3,3]); u=sym('U',[3,1]); Z = double(equationsToMatrix(equationsToMatrix(FC*u,FC(:)),u(:)));
            %     
            %     % evaluate numerically
            %     rng(1); norbits=4; ncenters=30; nsteps = 2; nvertices=2;
            %     u = rand(3,norbits,ncenters,nsteps); FC = rand(3,3,norbits); R_o2i = rand(3,3,norbits); 
            %     
            %     for i = 1:norbits;   iR_o2i(:,:,i) = inv(R_o2i(:,:,i)); end
            %     for i = 1:norbits; FC_proto(:,:,i) = iR_o2i(:,:,i)*FC(:,:,i)*R_o2i(:,:,i); end
            %     
            %     % explicit
            %     F1(1:3,1:ncenters*nsteps) = reshape(FC,3,3*norbits) * reshape(u,3*norbits,ncenters*nsteps);
            %     
            %     % using U matriz
            %     U = reshape(matmul_(Z,iR_o2i)   , 3, 9, 3, norbits, 1);
            %     U = matmul_(       reshape(R_o2i, 3, 3, 1, norbits, 1), U);
            %     U = sum_(     U .* reshape(u    , 1, 1, 3, norbits, ncenters*nsteps), 3);
            %     U =                 reshape(U   , 3, 3* 3* norbits, ncenters*nsteps);
            %     U = reshape(permute(U,[1,3,2])  , 3               * ncenters*nsteps, 9*norbits);
            %     
            %     X = U * FC_proto(:);
            %     
            %     U \ X - FC_proto(:)
            %     
            %     F1 - reshape(X,3,ncenters*nsteps)

            import am_lib.* am_dft.*

            %
            Z{1} = [];
            
            % the Z matrix factors out force constants leaving displacements
            % FC=sym('FC',[3,3]); U=sym('U',[3,1]); equationsToMatrix(equationsToMatrix(FC*U,FC(:)),U(:))
            Z{2} = [1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                    0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0;
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1].';
            
            % the Z matrix factors out force constants leaving displacements
            % THIS MATRIX IS PROBABLY WORNG. NEED TO DOUBLE CHECK.
            % THIS MATRIX IS PROBABLY WORNG. NEED TO DOUBLE CHECK.
            FC=sym('FC',[3,3,3]); U=sym('U',[3,3]); Z{3} = double(equationsToMatrix(equationsToMatrix(FC(:,:,1)*U+FC(:,:,2)*U+FC(:,:,3)*U,FC(:)),U(:)));
            % Matrix below generated with: FC=sym('FC',[3,3]); U=sym('U',[3,3]); Z{3} = double(equationsToMatrix(equationsToMatrix(FC*U,FC(:)),U(:)));
            Z{3} = [1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;
                    0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0;
                    0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1].';
              
            % get primitive pairs and link it to unit cell
            pp = get_primitive_cluster(ip,uc);
            
            % covert symmetries [pc-frac] -> [cart]
            sym_rebase_ = @(B,S) [[ matmul_(matmul_(B,S(1:3,1:3,:)),inv(B)), ...
                reshape(matmul_(B,S(1:3,4,:)),3,[],size(S,3))]; S(4,1:4,:)];
            Q_cart{1} = sym_rebase_(pp.bas,pp.Q{1});
            
            % initialize arrays
            nsteps = size(u,3);
            nFCs = sum(cellfun(@(x)size(x,2),ip.W));
            U = zeros(3*sum(nsteps*pp.ncenters),nFCs);
            I = zeros(3*sum(nsteps*pp.ncenters),1);
            X = reshape(1:numel(u),size(u));

            % record which cluster FCs belong to
            pc_natoms = numel(unique(pp.x2p));
            m_id = repelem([1:pc_natoms]   ,3*nsteps*pp.ncenters);
            s_id = repelem([1:ip.nclusters],cellfun(@(x)size(x,2),ip.W));

            % F [3m * 1] = - U [3m * nFcs] * FC [ nFCs * 1]: n pairs, m atoms ==> FC = - U \ F
            for m = 1:pc_natoms
                % get inds to properly reorder force constants
                I(m_id==m) = reshape(X(:,pp.c_id{m},:),[],1);
            for s = 1:ip.nclusters
                ex_ = [pp.ip_id{m}==s]; norbits= sum(ex_); nFCs = size(ip.W{s},2); ncenters = pp.ncenters(m);

                % define array reshape functions
               ushp_ = @(X) reshape(X, 3.^(ip.nvertices-1), norbits, ncenters*nsteps);
                shp_ = @(X) reshape(X, 3, 3.^ip.nvertices , norbits, ncenters*nsteps);
               fshp_ = @(X) reshape(X, 3,                      nFCs, ncenters*nsteps);

               % get rotation matrices
                R  = Q_cart{1}(1:3,1:3,pp.q_id{m}(ex_)); iR = permute(R,[2,1,3]);

                % setup displacement matrix
                switch ip.nvertices 
                    case 2
                        ux = u(:,pp.o_id{m}(ex_,:),:);
                    case 3
                        % construct displacement outer products
                        ux = outerc_( reshape(u(:,flatten_(pt.o{m}(ex_,:,1)),:),3,[]) , ...
                                      reshape(u(:,flatten_(pt.o{m}(ex_,:,2)),:),3,[]) );
                end
                
                % 0) reshape ux [3 x 3 x (npairs*natoms)] -> [9 x npairs x natoms]
                UW = ushp_(ux);
                % 1) transform outer dispacement product to irreducible orientation
                %    R^^2 [9 x 9] * UW [9 x 1  x npairs x natoms] = UW [9  x 1  x npairs x natoms]
                UW = matmul_(kronpow_(R,ip.nvertices-1),permute(UW,[1,4,2,3]));
                % 2) factor out force constants, leaving displacements
                %    Z [81 x 9] * UW [9 x 1  x npairs x natoms] =   UW [81 x 1  x npairs x natoms] reshaped to UW [3  x 27 x npairs x natoms]
                UW =  shp_(matmul_(Z{ip.nvertices},UW));
                % 3) return displacement to bond orientation
                %    iR [3 x 3] * UW [3 x 27 x npairs x natoms] =   UW [3  x 27 x npairs x natoms]
                UW =  shp_(matmul_(iR,UW));
                % 4) take into account intrinsic and crystallographic symmetries
                %    UW [3 x 27 x npairs x natoms] * bvt.W [ 27 * nfcs ] = UW [3 x ncfs x npairs x natoms] sum over pairs -> UW [3 x ncfs x 1 x natoms]
                UW = fshp_(matmul_(sum(UW,3),ip.W{s}));
                % 5) construct U
                U(m_id==m,s_id==s) = reshape(permute(UW,[1,3,2]), 3*ncenters*nsteps, nFCs );
                
            end
            end
        end

        function [q2u]   = expand_bvk_eigenvectors(bvk, uc, bz)
            % Expand primitive cell eigenvectors onto the unit cell; bz
            % must be the full mesh with the same dimensions as the
            % super unitcell (relative to that of the primitive cell).
            %
            % Use the expanded eigenvectors to convert normal phonon coordinates to
            % displacements and velocities (vectorized Wallace Eq. 10.41, p. 113):
            % U = q2u [ 1:3 * uc.natoms , bvk.nbranches * ibz.nks ]  * q_ks
            %
            % Note: Hermitian matrices eigendecomposed as A = U E U* for
            % a unitary matrix U and eigenvalues E -- meaning:
            %       q2u'*q2u = q2u*q2u' = Identity
            %
            % Q: Are the eigenvectors of a hermitian matrix unitary?
            % A: It appears to be so. Check with:
            %
            %        for i = 1:bz.nks;
            %               uni(i) = max(max(abs(inv(bz.U(:,:,i))-bz.U(:,:,i)')));
            %        end; semilogy(uni)
            %
            % Q: Does U diagonalize D?
            % A: Yes, U does diagonalize D when D is constructed properly including
            %    1/sqrt(M) factors! Otherwise it will leave it in block diagonal form
            %    with each block spanning the dimensions of the primitive cell basis.
            %    Check with:
            %
            %         n=[5;5;5];
            %         [bz,~] = get_zones(pc,n,''); bz = get_bvk_dispersion(bvk,bz);
            %         U   = expand_bvk_eigenvectors(bvk,uc,bz); N = normc_(real(U));
            %         q2u = real(U)  ./ N   ./ repelem(sqrt(uc.mass(uc.species)).',3,1);
            %         u2q = real(U)' ./ N.' .* repelem(sqrt(uc.mass(uc.species)).',3,1).';
            %
            %         M = repelem(uc.mass(uc.species),1,3);
            %         D = zeros(3*uc.natoms,3*uc.natoms);
            %         for m = 1:pp.pc_natoms
            %             for i = 1:pp.ncenters(m)
            %             for j = 1:pp.npairs(m)
            %                 x  = pp.c{m}(i);   xp = [1:3]+3*(x-1); ir = pp.i{m}(j);
            %                 y  = pp.o{m}(j,i); yp = [1:3]+3*(y-1); iq = pp.iq{m}(j);
            %
            %                 iphi = reshape(bvk.W{ir}*bvk.fc{ir}(:),3,3);
            %
            %                 D(xp,yp) = D(xp,yp) + ...
            %                     pp.Q{1}(1:3,1:3,iq) * permute(iphi,pp.Q{2}(:,iq)) * pp.Q{1}(1:3,1:3,iq)';
            %             end
            %             end
            %         end
            %         D = diag(1./sqrt(M)) * D * diag(1./sqrt(M));
            %         spy(abs(U'*D*U)>1E-8); view(2);  axis tight;
            %         spy(abs(u2q*D*q2u)>1E-4); view(2);  axis tight;
            %
            % Q: Are the eigenvalues obtained from the diagonalization the same?
            % A: The values are the same as that obtained by diagonalizing the
            %    primitive cell dynamical matrix. Eigenvenvalues obtained from the unit
            %    cell dynamical matrix are the same but sorted in a different order.
            %    Check with:
            %
            %       plot(bz.hw(:)./am_lib.units_eV,real(sqrt(real(diag(U'*D*U)))),'.');
            %
            % Q: Does real(U) span a complete basis?
            % A: No. While U does span a complete basis, the rank of real(U) is reduced and
            %    the null space of real(U) is increased by the same amount relative to
            %    that of U, indicating that real(U) vectors overlap on a smaller manifold.
            %

            import am_lib.* am_dft.*

            % orthonormalize U
            % U = bz.U;
            % for i = 1:bz.nks
            %     [~,~,a]=unique(rnd_(bz.hw(:,i)));
            %     for j = 1:max(a)
            %         U(:,a==j,i)=orth(U(:,a==j,i));
            %     end
            % end

            % find closest primitive lattice vector
            % G_ = @(tau) tau - inv(uc.tau2pc)*mod_(uc.tau2pc*tau);
            G_ = @(tau) tau;

            % expand primitive eigenvectors to unit cell
            % NOTE: exp in must ABSTOLUTELY be positive, i.e. +2i*pi and NOT -2i*pi
            u2p = flatten_([1:3].'+3*(uc.u2p-1));
            U = reshape(bz.U(u2p,:,:),3*uc.natoms,bvk.nbranches*bz.nks);
            E = repelem(exp(2i*pi*(uc.tau2pc*G_(uc.tau)).'*bz.k),3,bvk.nbranches);
            q2u = U .* E * sqrt(bvk.natoms/uc.natoms);

            % % explicit
            % q2u = zeros(uc.natoms,bvk.nbranches*bz.nks);
            % for i = 1:uc.natoms
            %     y = 0;  x = i; xp = [1:3] + 3*(x-1);
            %     m = uc.u2p(x); mp = [1:3] + 3*(m-1);
            % for j = 1:bz.nks; for k = 1:bvk.nbranches; y = y+1;
            %     q2u(xp,y) = bz.U(mp,k,j) * exp(-2i*pi*dot(uc.tau2pc*G_(uc.tau(:,x)),bz.k(:,j))) * sqrt(bvk.natoms/uc.natoms);
            % end;end
            % end
        end

        
        % aux phonons (anharmonic)

        function [f]   = get_bvt_forces(bvt,pt,u,algo)
            %
            % Get forces f [cart] from the displacement u [cart].
            %
            % Selecting a METHOD is optional (defaults to METHOD 2, which is ~ 10x faster):
            %
            % METHOD 1: Build U matrix and compute forces using:
            %           F [3m * 1] = - U [3m * nFcs] * FC [ nFCs * 1] for n pairs, m atoms
            %
            % METHOD 2: Compute F directly from tensor double dot products using:
            %           F [3 * m] = - FC [3 * 9n] * U [9n * m]:
            %           for n pairs, m atoms (9 ispans  outer displacement product)
            %
            % Compare forces obtained from the two methods:
            %     u = matmul_( md.bas, mod_( md.tau-uc.tau +.5 )-.5 );
            %     tic; f1 = get_bvt_forces(bvt,pt,u,1); toc
            %     tic; f2 = get_bvt_forces(bvt,pt,u,2); toc
            %     plot(f1(:),f2(:),'.')
            %
            import am_lib.* am_dft.*

            if nargin ~= 4; algo=2; end

            switch algo
                case 1
                    % get U matrix and indexing
                    [U,I] = get_U_matrix(bvt,pt,u);
                    % solve for forces f [3m * 1] = - U [3m * nFcs] * FC [ nFCs * 1]: n pairs, m atoms
                    f(I) = - U * [bvt.fc{:}].';
                    % rearrange forces
                    f = reshape(f,size(u));
                case 2
                    % get sizes
                    [~,natoms,nsteps] = size(u);
                    % build anharmonic force constant matrices for bonds
                    for m = 1:bvt.natoms
                        phi{m} = zeros(3,9*pt.npairs(m));
                    for j = 1:pt.npairs(m)
                        % get indicies
                        i = pt.i{m}(j); iq = pt.iq{m}(j);
                        % get irrep force constant indicies
                        iphi = reshape(bvt.W{i}*bvt.fc{i}(:),3,3,3);
                        % rotate force constants from irrep to orbit (NOTE: permuting iphi here
                        % does nothing because iphi already incorporates intrinsic symmetries.)
                        phi{m}(1:3,[1:9]+9*(j-1)) = reshape( kronpow_(pt.Q{1}(1:3,1:3,iq),3) * iphi(:) , [9,3] ).';
                    end
                    end
                    % compute forces on every atom at every step
                    f = zeros(3,natoms,nsteps);
                    for m = 1:bvt.natoms
                        % construct displacement outer products
                        ux = outerc_( reshape(u(:,flatten_(pt.o{m}(:,:,1)),:),3,[]) , ...
                                      reshape(u(:,flatten_(pt.o{m}(:,:,2)),:),3,[]) );
                        % reshape ux [3 x 3 x (npairs*natoms)] -> [9 x npairs x ncenters x nsteps]
                        ux = reshape(ux,9*pt.npairs(m),pt.ncenters(m),nsteps);
                        % compute forces
                        f(1:3,pt.c{m},:) = - matmul_(phi{m} , ux);
                    end
            end
        end

        % aux electrons


        function E        = eval_energies_(tb,x,k)
            % get hamiltonians
            nks = size(k,2); E = zeros(tb.nbands,nks); % recbas = inv(tb.bas).';
            for m = 1:nks
                % build input
                % input = num2cell([x,(recbas*k(:,m)).']);
                input = num2cell([x,k(:,m).']);
                % evaluate H
                E(:,m) = real(eig(tb.H(input{:})));
            end
            % sort values
            E = sort(E);
        end

    end
    
    % elastic constants
    
    
    methods (Static)

        function [C] = cij(C,flag)
            % [c11,c12,c44]
            switch flag
                case 'cubic'
                    [c11,c12,c44]=deal(C(1),C(2),C(3));
                    C = [c11 c12 c12 0.0 0.0 0.0; ...
                         c12 c11 c12 0.0 0.0 0.0; ...
                         c12 c12 c11 0.0 0.0 0.0; ...
                         0.0 0.0 0.0 c44 0.0 0.0; ...
                         0.0 0.0 0.0 0.0 c44 0.0; ...
                         0.0 0.0 0.0 0.0 0.0 c44 ];
                case ''
                otherwise; error('unknown symmetry type');     
            end
        end

        function [C] = cijkl2cij(CC)
            C = zeros(6,6);
            for im=1:3;for jm=1:3;for km=1:3;for lm=1:3
                if ( CC(im,jm,km,lm) ~= 0.0) 
                    [iv,jv]=ijkl2ij_local(im,jm,km,lm) ;
                    C(iv,jv) = CC(im,jm,km,lm);
                end
            end; end; end; end
            function [iv,jv] = ijkl2ij_local(ii,jj,kk,ll)   
                if (ii==1 && jj==1) iv=1; end
                if (ii==1 && jj==2) iv=6; end
                if (ii==1 && jj==3) iv=5; end
                if (ii==2 && jj==1) iv=6; end
                if (ii==2 && jj==2) iv=2; end
                if (ii==2 && jj==3) iv=4; end
                if (ii==3 && jj==1) iv=5; end
                if (ii==3 && jj==2) iv=4; end
                if (ii==3 && jj==3) iv=3; end
                if (kk==1 && ll==1) jv=1; end
                if (kk==1 && ll==2) jv=6; end
                if (kk==1 && ll==3) jv=5; end
                if (kk==2 && ll==1) jv=6; end
                if (kk==2 && ll==2) jv=2; end
                if (kk==2 && ll==3) jv=4; end
                if (kk==3 && ll==1) jv=5; end
                if (kk==3 && ll==2) jv=4; end
                if (kk==3 && ll==3) jv=3; end
            end
        end 
    end
    
    % unix functions and scripts

    methods (Static)

        % stand-alone

        function generate_scripts()

            import am_lib.* am_dft.*

            % initialize counter
            fprintf('Generating scripts:\n')

            % scripts functions
            fname='outcar2fp.sh'; fid=fopen([fname],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); fprintf(' ... %s (succeeded)\n',fname);
%{
#!/bin/bash
# Parse outcar to extract forces and displacements at each timestep.
# Antonio Mei Nov/2014
# Antonio Mei Jan/2017
usage_ () {
echo "Creates infile.force_position based on supplied outcar files."
echo ""
echo "Usage: $0 [-h] [-t] [-f] -o <outcar_list> -n <natoms> [-c <compress_name>]"
echo ""
echo "Example: $0 -f -t -o \"\$(find . -name "OUTCAR*" | grep 4.00-300)\" -n 250"
echo ""
echo "-h : prints this message"
echo "-n : [REQUIRED] number of atoms in the simulation cell"
echo "-o : [REQUIRED] list of outcar files to parse"
echo "-t : trims the last md run (useful for removing runs which have not completed)"
echo "-f : overwrites existing infile.force_position"
echo "-c : compresses infile.force_position to a tar.gz file"
echo ""
echo "infile.force_position file contents:"
echo "   x position   y position   z position     x force      y force      z force"
exit 1
}

main_ () {
# trim the last md run which may not have completed
trim_ () { tac $1 | awk '!found && /POSITION/{found=1;next}1' | tac ; }
# get position and forces
get_  () { cat $2 | grep -h -A $(($1+1)) POSITION  ; }
# cut header lines
cut_  () { cat $1 | sed '/^--$/d' | sed '/--/d' | sed '/POSITION/d' ; }
# compress produced infile.force_position
compress_ () { tar -zcvf infile.force_position.tar.gz infile.force_position ; }
#
if ${ISFORCE}; then
    if [ -f "./infile.force_position" ]; then
        rm ./infile.force_position
        printf " ... ./infile.force_position overwritten\n"
    fi
fi
#
if ${ISTRIM}; then
    printf " ... trim:\n"
    for F in "${FLIST}"; do
        printf " ...     %-100s\n" "${F}"
        trim_ ${F} | get_ ${NATOMS} | cut_ >> infile.force_position
    done
else
    printf " ... batch parsing without trim\n"
    get_ ${NATOMS} "${FLIST}" | cut_ >> infile.force_position
fi
#
printf " ... infile.force_position created\n"
#
if ${ISCOMPRESS}; then
    printf " ... infile.force_position.tar.gz compressed\n"
    compress_
fi
}

ISCOMPRESS=false; ISTRIM=false; ISFORCE=false;
if (($# == 0)); then usage_; exit 1; fi
while getopts "n:o:htfc" o; do
case "${o}" in
    o)  FLIST=${OPTARG} ;;
    n)  NATOMS=${OPTARG} ;;
    c)  ISCOMPRESS=true ;;
    t)  ISTRIM=true ;;
    f)  ISFORCE=true ;;
    h)  usage_; exit 0 ;;
    *)  usage_; exit 1 ;;
esac
done
main_
%}

            fname='outcar2en.sh'; fid=fopen([fname],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); fprintf(' ... %s (succeeded)\n',fname);
%{
#!/bin/bash
# Preprocess outcar to remove the last run, which may not have finished
# Antonio Mei Nov/2014
# Antonio Mei Jan/2017
usage_ () {
    echo "Creates infile.electron_energies based on supplied outcar files."
    echo ""
    echo "Usage: $0 [-h] [-t] [-f] -o <outcar_list> -n <nbands> [-c <compress_name>]"
    echo ""
    echo "Example: $0 -f -t -o \"\$(find . -name "OUTCAR*" | grep 4.00-300)\" -n 751"
    echo ""
    echo "-h : prints this message"
    echo "-n : [REQUIRED] number of bands"
    echo "-o : [REQUIRED] list of outcar files to parse"
    echo "-t : trims the last md run (useful for removing runs which have not completed)"
    echo "-f : overwrites existing infile.electron_energies"
    echo "-c : compresses infile.electron_energies to a tar.gz file"
    echo ""
    echo "infile.electron_energies file contents:"
    echo "   n index    En energy    fn occupation"
    exit 1
}
main_ () {
    # trim the last md run which may not have completed
    trim_ () { tac $1 | awk '!found && /POSITION/{found=1;next}1' | tac ; }
    # get energies
    get_  () { cat $2 | grep -h -A ${1} occupation  ; }
    # cut header lines
    cut_  () { cat $1 | sed '/^--$/d' | sed '/--/d' | sed '/occupation/d' ; }
    # compress produced infile.electron_energies
    compress_ () { tar -zcvf infile.electron_energies.tar.gz infile.electron_energies ; }
    #
    if ${ISFORCE}; then
        if [ -f "./infile.electron_energies" ]; then
            rm ./infile.electron_energies
            printf " ... ./infile.electron_energies overwritten\n"
        fi
    fi
    #
    if ${ISTRIM}; then
        printf " ... trim:\n"
        for F in "${FLIST}"; do
            printf " ...     %-100s\n" "${F}"
            trim_ ${F} | get_ ${NBANDS} | cut_ >> infile.electron_energies
        done
    else
        printf " ... batch parsing without trim\n"
        get_ ${NBANDS} "${FLIST}" | cut_ >> infile.electron_energies
    fi
    #
    awk '{ print $2 }' infile.electron_energies > infile.electron_energies.tmp && mv infile.electron_energies.tmp infile.electron_energies
    #
    printf " ... infile.electron_energies created\n"
    #
    if ${ISCOMPRESS}; then
        printf " ... infile.electron_energies.tar.gz compressed\n"
        compress_
    fi
}
ISCOMPRESS=false; ISTRIM=false; ISFORCE=false;
if (($# == 0)); then usage_; exit 1; fi
while getopts "n:o:htfc" o; do
    case "${o}" in
        o)  FLIST=${OPTARG} ;;
        n)  NBANDS=${OPTARG} ;;
        c)  ISCOMPRESS=true ;;
        t)  ISTRIM=true ;;
        f)  ISFORCE=true ;;
        h)  usage_; exit 0 ;;
        *)  usage_; exit 1 ;;
    esac
done
main_
%}

            fname='queue_vasp.sh'; fid=fopen([fname],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); fprintf(' ... %s (succeeded)\n',fname);
%{
#!/bin/bash
# Antonio Mei Sep/2017
usage_ () {
    echo "Queues unfinished calculations based on folders containing INCAR,POSCAR,KPOINTS,POTCAR,QSUB."
    echo ""
    echo "Usage: $0 [-h] [-l] [-v] [-n] [-t:target]" 
    echo ""
    echo ""
    echo "-h : prints this message"
    echo "-l : runs on local computer without queuing"
    echo "-t : target file to match"
    echo "-v : shows directories with an incomplete OUTCAR"
    echo "-n : (dry run) prints directors which would get queued"
    echo ""
    exit 1
}

main_ () {
    odir=$PWD
    dlist=$(find . -name ${TARGET} -exec sh -c '(cd `dirname {}` && pwd  )' ';' )
    for d in ${dlist}; do
    cd ${d}
     if ! grep -sq "Total CPU time" OUTCAR; then
      if [[ -e INCAR   ]]; then
       if [[ -e POSCAR  ]]; then
        if [[ -e KPOINTS ]]; then
         if [[ -e POTCAR  ]]; then
          echo $PWD
          if ! ${ISDRY}; then
              if hash qsub   2>/dev/null; then qsub   QSUB; fi
              if hash sbatch 2>/dev/null; then sbatch QSUB; fi
              if local; then vasp; fi
          fi
         fi
        fi
       fi
      fi
     fi
    cd ${odir}
    done
}

show_ () {
    grep -LR "Voluntary" --include="OUTCAR" . 
}

ISDRY=false; ISLOCAL=false; TARGET='QSUB'
if (($# == 0)); then usage_; exit 1; fi
while getopts "nlvht:" o; do
    case "${o}" in
    t)  TARGET=${OPTARG} ;;
        n)  ISDRY=true ;;
        l)  ISLOCAL=true ;;
        v)  show_; exit 0 ;;
        h)  usage_; exit 0 ;;
        *)  usage_; exit 1 ;;
    esac
done
main_
%}
          
            fname='vasp_raman.py'; fid=fopen([fname],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); fprintf(' ... %s (succeeded)\n',fname);
%{
#!/usr/bin/env python
#
# vasp_raman.py v. 0.6.0
#
# Raman off-resonant activity calculator
# using VASP as a back-end.
#
# Contributors: Alexandr Fonari (Georgia Tech)
# Shannon Stauffer (UT Austin)
#
# URL: http://raman-sc.github.io
#
# MIT license, 2013 - 2016
#


import re
import sys


def MAT_m_VEC(m, v):
p = [ 0.0 for i in range(len(v)) ]
for i in range(len(m)):
    assert len(v) == len(m[i]), 'Length of the matrix row is not equal to the length of the vector'
    p[i] = sum( [ m[i][j]*v[j] for j in range(len(v)) ] )
return p


def T(m):
p = [[ m[i][j] for i in range(len( m[j] )) ] for j in range(len( m )) ]
return p


def parse_poscar(poscar_fh):
# modified subroutine from phonopy 1.8.3 (New BSD license)
#
poscar_fh.seek(0) # just in case
lines = poscar_fh.readlines()
#
scale = float(lines[1])
if scale < 0.0:
    print "[parse_poscar]: ERROR negative scale not implemented."
    sys.exit(1)
#
b = []
for i in range(2, 5):
    b.append([float(x)*scale for x in lines[i].split()[:3]])
#
vol = b[0][0]*b[1][1]*b[2][2] + b[1][0]*b[2][1]*b[0][2] + b[2][0]*b[0][1]*b[1][2] - \
      b[0][2]*b[1][1]*b[2][0] - b[2][1]*b[1][2]*b[0][0] - b[2][2]*b[0][1]*b[1][0]
#
try:
    num_atoms = [int(x) for x in lines[5].split()]
    line_at = 6
except ValueError:
    symbols = [x for x in lines[5].split()]
    num_atoms = [int(x) for x in lines[6].split()]
    line_at = 7
nat = sum(num_atoms)
#
if lines[line_at][0].lower() == 's':
    line_at += 1
#
if (lines[line_at][0].lower() == 'c' or lines[line_at][0].lower() == 'k'):
    is_scaled = False
else:
    is_scaled = True
#
line_at += 1
#
positions = []
for i in range(line_at, line_at + nat):
    pos = [float(x) for x in lines[i].split()[:3]]
    #
    if is_scaled:
        pos = MAT_m_VEC(T(b), pos)
    #
    positions.append(pos)
#
poscar_header = ''.join(lines[1:line_at-1]) # will add title and 'Cartesian' later
return nat, vol, b, positions, poscar_header


def parse_env_params(params):
tmp = params.strip().split('_')
if len(tmp) != 4:
    print "[parse_env_params]: ERROR there should be exactly four parameters"
    sys.exit(1)
#
[first, last, nderiv, step_size] = [int(tmp[0]), int(tmp[1]), int(tmp[2]), float(tmp[3])]
#
return first, last, nderiv, step_size


#### subs for the output from VTST tools
def parse_freqdat(freqdat_fh, nat):
freqdat_fh.seek(0) # just in case
#
eigvals = [ 0.0 for i in range(nat*3) ]
#
for i in range(nat*3): # all frequencies should be supplied, regardless of requested to calculate
    tmp = freqdat_fh.readline().split()
    eigvals[i] = float(tmp[0])
#
return eigvals
#
def parse_modesdat(modesdat_fh, nat):
from math import sqrt
modesdat_fh.seek(0) # just in case
#
eigvecs = [ 0.0 for i in range(nat*3) ]
norms =   [ 0.0 for i in range(nat*3) ]
#
for i in range(nat*3): # all frequencies should be supplied, regardless of requested to calculate
    eigvec = []
    for j in range(nat):
        tmp = modesdat_fh.readline().split()
        eigvec.append([ float(tmp[x]) for x in range(3) ])
    #
    modesdat_fh.readline().split() # empty line
    eigvecs[i] = eigvec
    norms[i] = sqrt( sum( [abs(x)**2 for sublist in eigvec for x in sublist] ) )
#
return eigvecs, norms
#### end subs for VTST
#
def get_modes_from_OUTCAR(outcar_fh, nat):
from math import sqrt
eigvals = [ 0.0 for i in range(nat*3) ]
eigvecs = [ 0.0 for i in range(nat*3) ]
norms   = [ 0.0 for i in range(nat*3) ]
#
outcar_fh.seek(0) # just in case
while True:
    line = outcar_fh.readline()
    if not line:
        break
    #
    if "Eigenvectors after division by SQRT(mass)" in line:
        outcar_fh.readline() # empty line
        outcar_fh.readline() # Eigenvectors and eigenvalues of the dynamical matrix
        outcar_fh.readline() # ----------------------------------------------------
        outcar_fh.readline() # empty line
        #
        for i in range(nat*3): # all frequencies should be supplied, regardless of those requested to calculate
            outcar_fh.readline() # empty line
            p = re.search(r'^\s*(\d+).+?([\.\d]+) cm-1', outcar_fh.readline())
            eigvals[i] = float(p.group(2))
            #
            outcar_fh.readline() # X         Y         Z           dx          dy          dz
            eigvec = []
            #
            for j in range(nat):
                tmp = outcar_fh.readline().split()
                #
                eigvec.append([ float(tmp[x]) for x in range(3,6) ])
                #
            eigvecs[i] = eigvec
            norms[i] = sqrt( sum( [abs(x)**2 for sublist in eigvec for x in sublist] ) )
        #
        return eigvals, eigvecs, norms
    #
print "[get_modes_from_OUTCAR]: ERROR Couldn't find 'Eigenvectors after division by SQRT(mass)' in OUTCAR. Use 'NWRITE=3' in INCAR. Exiting..."
sys.exit(1)
#
def get_epsilon_from_OUTCAR(outcar_fh):
epsilon = []
#
outcar_fh.seek(0) # just in case
while True:
    line = outcar_fh.readline()
    if not line:
        break
    #
    if "HEAD OF MICROSCOPIC STATIC DIELECTRIC TENSOR (INDEPENDENT PARTICLE, excluding Hartree and local field effects)" in line:
        outcar_fh.readline()
        epsilon.append([float(x) for x in outcar_fh.readline().split()])
        epsilon.append([float(x) for x in outcar_fh.readline().split()])
        epsilon.append([float(x) for x in outcar_fh.readline().split()])
        return epsilon
#
raise RuntimeError("[get_epsilon_from_OUTCAR]: ERROR Couldn't find dielectric tensor in OUTCAR")
return 1
#
if __name__ == '__main__':
from math import pi
from shutil import move
import os
import datetime
import time
#import argparse
import optparse
#
print ""
print "    Raman off-resonant activity calculator,"
print "    using VASP as a back-end."
print ""
print "    Contributors: Alexandr Fonari  (Georgia Tech)"
print "                  Shannon Stauffer (UT Austin)"
print "    MIT License, 2013"
print "    URL: http://raman-sc.github.io"
print "    Started at: "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
print ""
#
description  = "Before run, set environment variables:\n"
description += "    VASP_RAMAN_RUN='mpirun vasp'\n"
description += "    VASP_RAMAN_PARAMS='[first-mode]_[last-mode]_[nderiv]_[step-size]'\n\n"
description += "bash one-liner is:\n"
description += "VASP_RAMAN_RUN='mpirun vasp' VASP_RAMAN_PARAMS='1_2_2_0.01' python vasp_raman.py"
#
parser = optparse.OptionParser(description=description)
parser.add_option('-g', '--gen', help='Generate POSCAR only', action='store_true')
parser.add_option('-u', '--use_poscar', help='Use provided POSCAR in the folder, USE WITH CAUTION!!', action='store_true')
(options, args) = parser.parse_args()
#args = vars(parser.parse_args())
args = vars(options)
#
VASP_RAMAN_RUN = os.environ.get('VASP_RAMAN_RUN')
if VASP_RAMAN_RUN == None:
    print "[__main__]: ERROR Set environment variable 'VASP_RAMAN_RUN'"
    print ""
    parser.print_help()
    sys.exit(1)
print "[__main__]: VASP_RAMAN_RUN='"+VASP_RAMAN_RUN+"'"
#
VASP_RAMAN_PARAMS = os.environ.get('VASP_RAMAN_PARAMS')
if VASP_RAMAN_PARAMS == None:
    print "[__main__]: ERROR Set environment variable 'VASP_RAMAN_PARAMS'"
    print ""
    parser.print_help()
    sys.exit(1)
print "[__main__]: VASP_RAMAN_PARAMS='"+VASP_RAMAN_PARAMS+"'"
#
first, last, nderiv, step_size = parse_env_params(VASP_RAMAN_PARAMS)
assert first >= 1,    '[__main__]: First mode should be equal or larger than 1'
assert last >= first, '[__main__]: Last mode should be equal or larger than first mode'
if args['gen']: assert last == first, "[__main__]: '-gen' mode -> only generation for the one mode makes sense"
assert nderiv == 2,   '[__main__]: At this time, nderiv = 2 is the only supported'
disps = [-1, 1]      # hardcoded for
coeffs = [-0.5, 0.5] # three point stencil (nderiv=2)
#
try:
    poscar_fh = open('POSCAR.phon', 'r')
except IOError:
    print "[__main__]: ERROR Couldn't open input file POSCAR.phon, exiting...\n"
    sys.exit(1)
#
# nat, vol, b, poscar_header = parse_poscar_header(poscar_fh)
nat, vol, b, pos, poscar_header = parse_poscar(poscar_fh)
print pos
#print poscar_header
#sys.exit(0)
#
# either use modes from vtst tools or VASP
if os.path.isfile('freq.dat') and os.path.isfile('modes_sqrt_amu.dat'):
    try:
        freqdat_fh = open('freq.dat', 'r')
    except IOError:
        print "[__main__]: ERROR Couldn't open freq.dat, exiting...\n"
        sys.exit(1)
    #
    eigvals = parse_freqdat(freqdat_fh, nat)
    freqdat_fh.close()
    #
    try: 
        modes_fh = open('modes_sqrt_amu.dat' , 'r')
    except IOError:
        print "[__main__]: ERROR Couldn't open modes_sqrt_amu.dat, exiting...\n"
        sys.exit(1)
    #
    eigvecs, norms = parse_modesdat(modes_fh, nat)
    modes_fh.close()
#
elif os.path.isfile('OUTCAR.phon'):
    try:
        outcar_fh = open('OUTCAR.phon', 'r')
    except IOError:
        print "[__main__]: ERROR Couldn't open OUTCAR.phon, exiting...\n"
        sys.exit(1)
    #
    eigvals, eigvecs, norms = get_modes_from_OUTCAR(outcar_fh, nat)
    outcar_fh.close()
#
else:
    print "[__main__]: Neither OUTCAR.phon nor freq.dat/modes_sqrt_amu.dat were found, nothing to do, exiting..."
    sys.exit(1)
#
output_fh = open('vasp_raman.dat', 'w')
output_fh.write("# mode    freq(cm-1)    alpha    beta2    activity\n")
for i in range(first-1, last):
    eigval = eigvals[i]
    eigvec = eigvecs[i]
    norm = norms[i]
    #
    print ""
    print "[__main__]: Mode #%i: frequency %10.7f cm-1; norm: %10.7f" % ( i+1, eigval, norm )
    #
    ra = [[0.0 for x in range(3)] for y in range(3)]
    for j in range(len(disps)):
        disp_filename = 'OUTCAR.%04d.%+d.out' % (i+1, disps[j])
        #
        try:
            outcar_fh = open(disp_filename, 'r')
            print "[__main__]: File "+disp_filename+" exists, parsing..."
        except IOError:
            if args['use_poscar'] != True:
                print "[__main__]: File "+disp_filename+" not found, preparing displaced POSCAR"
                poscar_fh = open('POSCAR', 'w')
                poscar_fh.write("%s %4.1e \n" % (disp_filename, step_size))
                poscar_fh.write(poscar_header)
                poscar_fh.write("Cartesian\n")
                #
                for k in range(nat):
                    pos_disp = [ pos[k][l] + eigvec[k][l]*step_size*disps[j]/norm for l in range(3)]
                    poscar_fh.write( '%15.10f %15.10f %15.10f\n' % (pos_disp[0], pos_disp[1], pos_disp[2]) )
                    #print '%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f' % (pos[k][0], pos[k][1], pos[k][2], dis[k][0], dis[k][1], dis[k][2])
                poscar_fh.close()
            else:
                print "[__main__]: Using provided POSCAR"
            #
            if args['gen']: # only generate POSCARs
                poscar_fn = 'POSCAR.%+d.out' % disps[j]
                move('POSCAR', poscar_fn)
                print "[__main__]: '-gen' mode -> "+poscar_fn+" with displaced atoms have been generated"
                #
                if j+1 == len(disps): # last iteration for the current displacements list
                    print "[__main__]: '-gen' mode -> POSCAR files with displaced atoms have been generated, exiting now"
                    sys.exit(0)
            else: # run VASP here
                print "[__main__]: Running VASP..."
                os.system(VASP_RAMAN_RUN)
                try:
                    move('OUTCAR', disp_filename)
                except IOError:
                    print "[__main__]: ERROR Couldn't find OUTCAR file, exiting..."
                    sys.exit(1)
                #
                outcar_fh = open(disp_filename, 'r')
        #
        try:
            eps = get_epsilon_from_OUTCAR(outcar_fh)
            outcar_fh.close()
        except Exception, err:
            print err
            print "[__main__]: Moving "+disp_filename+" back to 'OUTCAR' and exiting..."
            move(disp_filename, 'OUTCAR')
            sys.exit(1)
        #
        for m in range(3):
            for n in range(3):
                ra[m][n]   += eps[m][n] * coeffs[j]/step_size * norm * vol/(4.0*pi)
        #units: A^2/amu^1/2 =         dimless   * 1/A         * 1/amu^1/2  * A^3
    #
    alpha = (ra[0][0] + ra[1][1] + ra[2][2])/3.0
    beta2 = ( (ra[0][0] - ra[1][1])**2 + (ra[0][0] - ra[2][2])**2 + (ra[1][1] - ra[2][2])**2 + 6.0 * (ra[0][1]**2 + ra[0][2]**2 + ra[1][2]**2) )/2.0
    print ""
    print "! %4i  freq: %10.5f  alpha: %10.7f  beta2: %10.7f  activity: %10.7f " % (i+1, eigval, alpha, beta2, 45.0*alpha**2 + 7.0*beta2)
    output_fh.write("%03i  %10.5f  %10.7f  %10.7f  %10.7f\n" % (i+1, eigval, alpha, beta2, 45.0*alpha**2 + 7.0*beta2))
    output_fh.flush()
#
output_fh.close()
%}

        end

    end


    % mex functions

    methods (Static)

        function compile_mex()
            % compile mex functions
            %
            % to compile the library manually
            %
            %       fort -O3 -parallel -fpp -c am_mex_lib.f90
            %
            % to compile individual interfaces manually
            %
            %       ifort -O3 -parallel -fpp -fPIC -L/Applications/MATLAB_R2016b.app/bin/maci64
            %       -I/Applications/MATLAB_R2016b.app/extern/include -lmx -lmex -lmat -nofor_main
            %       -bundle ./uc2ws_mex.f90 -o uc2ws_mex.mexmaci64
            %

            import am_lib.* am_dft.*

            % define mex parameters
            FC        = 'ifort';
            FFLAGS    = '-O3 -parallel -fpp -fPIC -lmx -lmex -lmat -nofor_main -bundle -implicitnone -assume realloc_lhs';
            MPATH     = '/Applications/MATLAB_R2017b.app';
            LIBS      = ['-L',MPATH,'/bin/maci64 -I',MPATH,'/extern/include'];
            EXT       = '.mexmaci64';
            DEBUG     = '-debug';

            % initialize counter
            i=0;

            % mex library
            flib = 'am_mex_lib'; fid=fopen([flib,'.f90'],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid);

% START OF FORTRAN CODE
%{
module am_mex_lib

    implicit none

    public

    !  $$$$$$\   $$$$$$\  $$\   $$\  $$$$$$\ $$$$$$$$\  $$$$$$\  $$\   $$\ $$$$$$$$\  $$$$$$\
    ! $$  __$$\ $$  __$$\ $$$\  $$ |$$  __$$\\__$$  __|$$  __$$\ $$$\  $$ |\__$$  __|$$  __$$\
    ! $$ /  \__|$$ /  $$ |$$$$\ $$ |$$ /  \__|  $$ |   $$ /  $$ |$$$$\ $$ |   $$ |   $$ /  \__|
    ! $$ |      $$ |  $$ |$$ $$\$$ |\$$$$$$\    $$ |   $$$$$$$$ |$$ $$\$$ |   $$ |   \$$$$$$\
    ! $$ |      $$ |  $$ |$$ \$$$$ | \____$$\   $$ |   $$  __$$ |$$ \$$$$ |   $$ |    \____$$\
    ! $$ |  $$\ $$ |  $$ |$$ |\$$$ |$$\   $$ |  $$ |   $$ |  $$ |$$ |\$$$ |   $$ |   $$\   $$ |
    ! \$$$$$$  | $$$$$$  |$$ | \$$ |\$$$$$$  |  $$ |   $$ |  $$ |$$ | \$$ |   $$ |   \$$$$$$  |
    !  \______/  \______/ \__|  \__| \______/   \__|   \__|  \__|\__|  \__|   \__|    \______/

    integer    , parameter :: sp      = kind(0.0E0) !> single precision
    integer    , parameter :: dp      = kind(0.0D0) !> double precision
    real(dp)   , parameter :: tiny    =  1.0D-5
    real(dp)   , parameter :: halfpi  =  1.570796326794897_dp
    real(dp)   , parameter :: sqrtpi  =  1.772453850905516_dp
    real(dp)   , parameter :: pi      =  3.141592653589793_dp
    real(dp)   , parameter :: twopi   =  6.283185307179586_dp
    real(dp)   , parameter :: fourpi  = 12.566370614359172_dp
    complex(dp), parameter :: cmplx_i = cmplx(0,1,dp)
    complex(dp), parameter :: itwopi  = cmplx_i*twopi

    ! $$$$$$\ $$\   $$\ $$$$$$$$\ $$$$$$$$\ $$$$$$$\  $$$$$$$$\  $$$$$$\   $$$$$$\  $$$$$$$$\  $$$$$$\
    ! \_$$  _|$$$\  $$ |\__$$  __|$$  _____|$$  __$$\ $$  _____|$$  __$$\ $$  __$$\ $$  _____|$$  __$$\
    !   $$ |  $$$$\ $$ |   $$ |   $$ |      $$ |  $$ |$$ |      $$ /  $$ |$$ /  \__|$$ |      $$ /  \__|
    !   $$ |  $$ $$\$$ |   $$ |   $$$$$\    $$$$$$$  |$$$$$\    $$$$$$$$ |$$ |      $$$$$\    \$$$$$$\
    !   $$ |  $$ \$$$$ |   $$ |   $$  __|   $$  __$$< $$  __|   $$  __$$ |$$ |      $$  __|    \____$$\
    !   $$ |  $$ |\$$$ |   $$ |   $$ |      $$ |  $$ |$$ |      $$ |  $$ |$$ |  $$\ $$ |      $$\   $$ |
    ! $$$$$$\ $$ | \$$ |   $$ |   $$$$$$$$\ $$ |  $$ |$$ |      $$ |  $$ |\$$$$$$  |$$$$$$$$\ \$$$$$$  |
    ! \______|\__|  \__|   \__|   \________|\__|  \__|\__|      \__|  \__| \______/ \________| \______/

    interface foursort
        module procedure d_foursort, i_foursort
    end interface ! foursort

    interface fourrank
        module procedure d_fourrank, i_fourrank
    end interface ! fourrank

    contains

    !  $$$$$$\ $$$$$$$$\  $$$$$$\ $$$$$$$$\ $$$$$$\  $$$$$$\ $$$$$$$$\ $$$$$$\  $$$$$$\   $$$$$$\
    ! $$  __$$\\__$$  __|$$  __$$\\__$$  __|\_$$  _|$$  __$$\\__$$  __|\_$$  _|$$  __$$\ $$  __$$\
    ! $$ /  \__|  $$ |   $$ /  $$ |  $$ |     $$ |  $$ /  \__|  $$ |     $$ |  $$ /  \__|$$ /  \__|
    ! \$$$$$$\    $$ |   $$$$$$$$ |  $$ |     $$ |  \$$$$$$\    $$ |     $$ |  $$ |      \$$$$$$\
    !  \____$$\   $$ |   $$  __$$ |  $$ |     $$ |   \____$$\   $$ |     $$ |  $$ |       \____$$\
    ! $$\   $$ |  $$ |   $$ |  $$ |  $$ |     $$ |  $$\   $$ |  $$ |     $$ |  $$ |  $$\ $$\   $$ |
    ! \$$$$$$  |  $$ |   $$ |  $$ |  $$ |   $$$$$$\ \$$$$$$  |  $$ |   $$$$$$\ \$$$$$$  |\$$$$$$  |
    !  \______/   \__|   \__|  \__|  \__|   \______| \______/   \__|   \______| \______/  \______/

    pure function factorial(n) result(y)
        !
        implicit none
        !
        integer, intent(in) :: n
        real(dp) :: y
        !
        if (n.ge.1) then
            y = product([1:n])
        else
            y = 1
        endif
    end function  factorial

    pure function nchoosek(n,k) result (res)
        !
        implicit none
        !
        integer, intent(in) :: n
        integer, intent(in) :: k
        integer :: res
        !
        res=factorial(n)/(factorial(k)*factorial(n-k))
        !
    end function  nchoosek

    function      perms(n) result(PT)
        !
        implicit none
        !
        integer, intent(in) :: n
        integer, allocatable :: P(:,:)
        integer, allocatable :: PT(:,:)
        integer :: m
        !
        m = product([1:n])
        !
        allocate(P(m,n))
        !
        call permutate([1:n],P)
        !
        allocate(PT,source=transpose(P))
        !
        contains
        recursive subroutine permutate(E, P)
            !
            implicit none
            !
            integer, intent(in)  :: E(:) ! array of objects
            integer, intent(out) :: P(:,:) ! permutations of E
            integer :: N, Nfac, i, k, S(size(P,1)/size(E), size(E)-1)
            N = size(E); Nfac = size(P,1);
            do i = 1, N
              if( N>1 ) call permutate((/E(:i-1), E(i+1:)/), S)
              forall(k=1:Nfac/N) P((i-1)*Nfac/N+k,:) = (/E(i), S(k,:)/)
            enddo
        end subroutine permutate
    end function  perms

    ! $$\   $$\ $$\   $$\ $$\      $$\ $$$$$$$\  $$$$$$$$\ $$$$$$$\        $$$$$$$$\ $$\   $$\ $$$$$$$$\  $$$$$$\  $$$$$$$\ $$\     $$\
    ! $$$\  $$ |$$ |  $$ |$$$\    $$$ |$$  __$$\ $$  _____|$$  __$$\       \__$$  __|$$ |  $$ |$$  _____|$$  __$$\ $$  __$$\\$$\   $$  |
    ! $$$$\ $$ |$$ |  $$ |$$$$\  $$$$ |$$ |  $$ |$$ |      $$ |  $$ |         $$ |   $$ |  $$ |$$ |      $$ /  $$ |$$ |  $$ |\$$\ $$  /
    ! $$ $$\$$ |$$ |  $$ |$$\$$\$$ $$ |$$$$$$$\ |$$$$$\    $$$$$$$  |         $$ |   $$$$$$$$ |$$$$$\    $$ |  $$ |$$$$$$$  | \$$$$  /
    ! $$ \$$$$ |$$ |  $$ |$$ \$$$  $$ |$$  __$$\ $$  __|   $$  __$$<          $$ |   $$  __$$ |$$  __|   $$ |  $$ |$$  __$$<   \$$  /
    ! $$ |\$$$ |$$ |  $$ |$$ |\$  /$$ |$$ |  $$ |$$ |      $$ |  $$ |         $$ |   $$ |  $$ |$$ |      $$ |  $$ |$$ |  $$ |   $$ |
    ! $$ | \$$ |\$$$$$$  |$$ | \_/ $$ |$$$$$$$  |$$$$$$$$\ $$ |  $$ |         $$ |   $$ |  $$ |$$$$$$$$\  $$$$$$  |$$ |  $$ |   $$ |
    ! \__|  \__| \______/ \__|     \__|\_______/ \________|\__|  \__|         \__|   \__|  \__|\________| \______/ \__|  \__|   \__|

    pure function primes(nprimes)
        !
        ! naive approach to generating the first n prime numbers
        !
        implicit none
        !
        integer, intent(in)  :: nprimes
        integer, allocatable :: primes(:) ! array that will hold the primes
        integer :: at, found, i
        logical :: is_prime
        !
        allocate (primes(nprimes))
        !
        primes(1) = 2
        at = 2
        found = 1
        do
            is_prime = .true. ! assume prime
            do i = 1, found
                if (modulo(at,primes(i)).eq.0) then ! if divisible by any other element
                    is_prime = .false.               ! in the array, then not prime.
                    at = at + 1
                    continue
                end if
            end do
            found = found + 1
            primes(found) = at
            at = at + 1
            if (found == nprimes) then ! stop when all primes are found
                exit
            endif
        end do
        !
    end function  primes

    !  $$$$$$\   $$$$$$\  $$$$$$$\ $$$$$$$$\ $$$$$$\ $$\   $$\  $$$$$$\
    ! $$  __$$\ $$  __$$\ $$  __$$\\__$$  __|\_$$  _|$$$\  $$ |$$  __$$\
    ! $$ /  \__|$$ /  $$ |$$ |  $$ |  $$ |     $$ |  $$$$\ $$ |$$ /  \__|
    ! \$$$$$$\  $$ |  $$ |$$$$$$$  |  $$ |     $$ |  $$ $$\$$ |$$ |$$$$\
    !  \____$$\ $$ |  $$ |$$  __$$<   $$ |     $$ |  $$ \$$$$ |$$ |\_$$ |
    ! $$\   $$ |$$ |  $$ |$$ |  $$ |  $$ |     $$ |  $$ |\$$$ |$$ |  $$ |
    ! \$$$$$$  | $$$$$$  |$$ |  $$ |  $$ |   $$$$$$\ $$ | \$$ |\$$$$$$  |
    !  \______/  \______/ \__|  \__|  \__|   \______|\__|  \__| \______/


    pure function d_foursort(x) result(y)
        !
        implicit none
        !
        real(dp), intent(in) :: x(4)
        real(dp) :: y(4)
        !
        y = x
        if (y(1).gt.y(2)) then; y([1,2]) = y([2,1]); endif
        if (y(3).gt.y(4)) then; y([3,4]) = y([4,3]); endif
        if (y(1).gt.y(3)) then; y([1,3]) = y([3,1]); endif
        if (y(2).gt.y(4)) then; y([2,4]) = y([4,2]); endif
        if (y(2).gt.y(3)) then; y([2,3]) = y([3,2]); endif
        !
    end function  d_foursort

    pure function i_foursort(x) result(y)
        !
        implicit none
        !
        integer, intent(in) :: x(4)
        integer :: y(4)
        !
        y = x
        if (y(1).gt.y(2)) then; y([1,2]) = y([2,1]); endif
        if (y(3).gt.y(4)) then; y([3,4]) = y([4,3]); endif
        if (y(1).gt.y(3)) then; y([1,3]) = y([3,1]); endif
        if (y(2).gt.y(4)) then; y([2,4]) = y([4,2]); endif
        if (y(2).gt.y(3)) then; y([2,3]) = y([3,2]); endif
        !
    end function  i_foursort

    pure function d_fourrank(n) result(ind)
        ! borrowed from olle, who took it from stack overflow
        implicit none
        !
        real(dp),intent(in) :: n(4)
        integer :: ind(4)
        integer :: low1,high1,low2,high2,highest,lowest,middle1,middle2
        !
        if ( n(1) <= n(2) ) then
            low1 = 1
            high1 = 2
        else
            low1 = 2
            high1 = 1
        endif

        if ( n(3) <= n(4) ) then
            low2 = 3
            high2 = 4
        else
            low2 = 4
            high2 = 3
        endif

        if ( n(low1) <= n(low2) ) then
            lowest = low1
            middle1 = low2
        else
            lowest = low2
            middle1 = low1
        endif

        if ( n(high1) >= n(high2) ) then
            highest = high1
            middle2 = high2
        else
            highest = high2
            middle2 = high1
        endif

        if ( n(middle1) < n(middle2) ) then
            ind=(/lowest,middle1,middle2,highest/)
        else
            ind=(/lowest,middle2,middle1,highest/)
        endif
        !
    end function  d_fourrank

    pure function i_fourrank(n) result(ind)
        ! borrowed from olle, who took it from stack overflow
        implicit none
        !
        integer,intent(in) :: n(4)
        integer :: ind(4)
        integer :: low1,high1,low2,high2,highest,lowest,middle1,middle2
        !
        if ( n(1) <= n(2) ) then
            low1 = 1
            high1 = 2
        else
            low1 = 2
            high1 = 1
        endif

        if ( n(3) <= n(4) ) then
            low2 = 3
            high2 = 4
        else
            low2 = 4
            high2 = 3
        endif

        if ( n(low1) <= n(low2) ) then
            lowest = low1
            middle1 = low2
        else
            lowest = low2
            middle1 = low1
        endif

        if ( n(high1) >= n(high2) ) then
            highest = high1
            middle2 = high2
        else
            highest = high2
            middle2 = high1
        endif

        if ( n(middle1) < n(middle2) ) then
            ind=(/lowest,middle1,middle2,highest/)
        else
            ind=(/lowest,middle2,middle1,highest/)
        endif
        !
    end function  i_fourrank

    pure function d_rank(v) result(ind)
      !
      ! based on dlasrt2
      !
      implicit none
      !
      real(dp), intent(in) :: v(:)
      real(dp), allocatable :: d(:)
      integer, allocatable :: ind(:)
      integer :: select, endd, i, j, start, stkpnt, tmpind, stack(2,32), n
      real(dp) :: d1, d2, d3, dmnmx, tmp
      !
      select = 20
      !
      n = size(v)
      !
      allocate(d,source=v)
      !
      allocate(ind,source=[1:n])
      !
      if (n.le.1) return
      !
      stkpnt = 1
      stack(1,1) = 1
      stack(2,1) = n
      !
      do
         start = stack(1, stkpnt)
         endd = stack(2, stkpnt)
         stkpnt = stkpnt - 1
         if(endd-start.gt.0) then
            do i = start + 1, endd
               do j = i, start + 1, -1
                  if(d(j).lt.d(j-1)) then
                     dmnmx = d(j)
                     d(j) = d(j-1)
                     d(j-1) = dmnmx
                     tmpind = ind(j)
                     ind(j) = ind(j-1)
                     ind(j-1) = tmpind
                  else
                     exit
                  end if
               enddo
            enddo
         else if(endd-start.gt.select) then
            d1 = d(start)
            d2 = d(endd)
            i = (start+endd) / 2
            d3 = d(i)
            if(d1.lt.d2) then
               if(d3.lt.d1) then
                  dmnmx = d1
               else if(d3.lt.d2) then
                  dmnmx = d3
               else
                  dmnmx = d2
               end if
            else
               if(d3.lt.d2) then
                  dmnmx = d2
               else if(d3.lt.d1) then
                  dmnmx = d3
               else
                  dmnmx = d1
               end if
            end if
            i = start - 1
            j = endd + 1
            do
               do
                  j = j - 1
                  if (.not.(d(j).gt.dmnmx)) exit
               enddo
               do
                  i = i + 1
                  if (.not.(d(i).lt.dmnmx)) exit
               enddo
               if(i.gt.j) then
                  exit
               else
                  tmp = d(i)
                  d(i) = d(j)
                  d(j) = tmp
                  tmpind = ind(j)
                  ind(j) = ind(i)
                  ind(i) = tmpind
               end if
            enddo
            if(j-start.gt.endd-j-1) then
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = start
               stack(2, stkpnt) = j
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = j + 1
               stack(2, stkpnt) = endd
            else
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = j + 1
               stack(2, stkpnt) = endd
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = start
               stack(2, stkpnt) = j
            end if
         end if
         if (.not.(stkpnt.gt.0)) exit
      enddo
    end function  d_rank

    pure function i_rank(v) result(ind)
      !
      ! based on dlasrt2
      !
      implicit none
      !
      integer, intent(in) :: v(:)
      integer, allocatable :: d(:)
      integer, allocatable :: ind(:)
      integer :: select, endd, i, j, start, stkpnt, tmpind, stack(2,32), n
      integer :: d1, d2, d3, dmnmx, tmp
      !
      select = 20
      !
      n = size(v)
      !
      allocate(d,source=v)
      !
      allocate(ind,source=[1:n])
      !
      if (n.le.1) return
      !
      stkpnt = 1
      stack(1,1) = 1
      stack(2,1) = n
      !
      do
         start = stack(1, stkpnt)
         endd = stack(2, stkpnt)
         stkpnt = stkpnt - 1
         if(endd-start.gt.0) then
            do i = start + 1, endd
               do j = i, start + 1, -1
                  if(d(j).lt.d(j-1)) then
                     dmnmx = d(j)
                     d(j) = d(j-1)
                     d(j-1) = dmnmx
                     tmpind = ind(j)
                     ind(j) = ind(j-1)
                     ind(j-1) = tmpind
                  else
                     exit
                  end if
               enddo
            enddo
         else if(endd-start.gt.select) then
            d1 = d(start)
            d2 = d(endd)
            i = (start+endd) / 2
            d3 = d(i)
            if(d1.lt.d2) then
               if(d3.lt.d1) then
                  dmnmx = d1
               else if(d3.lt.d2) then
                  dmnmx = d3
               else
                  dmnmx = d2
               end if
            else
               if(d3.lt.d2) then
                  dmnmx = d2
               else if(d3.lt.d1) then
                  dmnmx = d3
               else
                  dmnmx = d1
               end if
            end if
            i = start - 1
            j = endd + 1
            do
               do
                  j = j - 1
                  if (.not.(d(j).gt.dmnmx)) exit
               enddo
               do
                  i = i + 1
                  if (.not.(d(i).lt.dmnmx)) exit
               enddo
               if(i.gt.j) then
                  exit
               else
                  tmp = d(i)
                  d(i) = d(j)
                  d(j) = tmp
                  tmpind = ind(j)
                  ind(j) = ind(i)
                  ind(i) = tmpind
               end if
            enddo
            if(j-start.gt.endd-j-1) then
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = start
               stack(2, stkpnt) = j
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = j + 1
               stack(2, stkpnt) = endd
            else
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = j + 1
               stack(2, stkpnt) = endd
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = start
               stack(2, stkpnt) = j
            end if
         end if
         if (.not.(stkpnt.gt.0)) exit
      enddo
    end function  i_rank

    ! $$\      $$\  $$$$$$\ $$$$$$$$\ $$\   $$\       $$$$$$$$\ $$\   $$\ $$\   $$\  $$$$$$\ $$$$$$$$\ $$$$$$\  $$$$$$\  $$\   $$\  $$$$$$\
    ! $$$\    $$$ |$$  __$$\\__$$  __|$$ |  $$ |      $$  _____|$$ |  $$ |$$$\  $$ |$$  __$$\\__$$  __|\_$$  _|$$  __$$\ $$$\  $$ |$$  __$$\
    ! $$$$\  $$$$ |$$ /  $$ |  $$ |   $$ |  $$ |      $$ |      $$ |  $$ |$$$$\ $$ |$$ /  \__|  $$ |     $$ |  $$ /  $$ |$$$$\ $$ |$$ /  \__|
    ! $$\$$\$$ $$ |$$$$$$$$ |  $$ |   $$$$$$$$ |      $$$$$\    $$ |  $$ |$$ $$\$$ |$$ |        $$ |     $$ |  $$ |  $$ |$$ $$\$$ |\$$$$$$\
    ! $$ \$$$  $$ |$$  __$$ |  $$ |   $$  __$$ |      $$  __|   $$ |  $$ |$$ \$$$$ |$$ |        $$ |     $$ |  $$ |  $$ |$$ \$$$$ | \____$$\
    ! $$ |\$  /$$ |$$ |  $$ |  $$ |   $$ |  $$ |      $$ |      $$ |  $$ |$$ |\$$$ |$$ |  $$\   $$ |     $$ |  $$ |  $$ |$$ |\$$$ |$$\   $$ |
    ! $$ | \_/ $$ |$$ |  $$ |  $$ |   $$ |  $$ |      $$ |      \$$$$$$  |$$ | \$$ |\$$$$$$  |  $$ |   $$$$$$\  $$$$$$  |$$ | \$$ |\$$$$$$  |
    ! \__|     \__|\__|  \__|  \__|   \__|  \__|      \__|       \______/ \__|  \__| \______/   \__|   \______| \______/ \__|  \__| \______/

    ! special functions

    pure function legendre(l,m,x) result(y)
        !
        ! Associated legrende polynomial
        !
        ! W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, Numerical Recipes in Fortran 77: The Art
        ! of Scientific Computing, 2 edition (Cambridge University Press, Cambridge England?; New York, 1992), p 246.
        !
        implicit none
        !
        integer , intent(in) :: l,m
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: y(:)
        integer  :: ll
        real(dp) :: pll(size(x))
        real(dp) :: pmm(size(x))
        real(dp) :: pmmp1(size(x))
        real(dp) :: somx2(size(x))
        !
        allocate(y(size(x)))
        !
        if ((m.lt.0).or.(m.gt.l).or.(all(abs(x).gt.1.0_dp))) then
            y = 0
        endif
        !
        pmm=1.0
        if (m > 0) then
            somx2=sqrt((1.0_dp-x)*(1.0_dp+x))
            pmm=product(arth(1.0_dp,2.0_dp,m))*somx2**m
            if (mod(m,2) == 1) pmm=-pmm
        end if
        if (l == m) then
            y=pmm
        else
            pmmp1=x*(2*m+1)*pmm
            if (l == m+1) then
                y=pmmp1
            else
                do ll=m+2,l
                    pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
                    pmm=pmmp1
                    pmmp1=pll
                end do
                y=pll
            end if
        end if
        contains
        pure function arth(first,increment,n)
            !
            implicit none
            !
            real(dp), intent(in) :: first
            real(dp), intent(in) :: increment
            integer , intent(in) :: n
            real(dp) :: arth(n)
            integer  :: k,k2
            real(dp) :: temp
            !
            if (n > 0) arth(1)=first
            if (n <= 16) then
                do k=2,n
                    arth(k)=arth(k-1)+increment
                enddo
            else
                do k=2,8
                    arth(k)=arth(k-1)+increment
                end do
                temp=increment*8
                k=8
                do
                    if (k >= n) exit
                    k2=k+k
                    arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
                    temp=temp+temp
                    k=k2
                enddo
            end if
        end function arth
    end function  legendre

    pure function laguerre(k,p,x) result(y)
        ! associated Laguerre polynomial L_k^p(x)(k,p,x)
        ! Using the expression from Samuel Shaw Ming Wong "Computational Methods in Physics and Engineering", Eq 4-86 p 139
        ! Note, there is a typographical error on http://mathworld.wolfram.com/AssociatedLaguerrePolynomial.html Eq 10
        ! Also see Linus Pauling "Introuction to Quantum Mechancs"
        implicit none
        !
        integer , intent(in) :: k, p
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: y(:)
        integer :: j
        !
        allocate(y,mold=x)
        y=0.0_dp
        !
        do j = 0, k
            y = y + nchoosek(k+p,k-j) * (-x)**j / factorial(j)
        enddo
        !
    end function  laguerre

    pure function Ylm(l,m,theta,phi)
        !
        ! computes spherical harmonics. Theta and phi in radians. size(theta) must equal size(phi); together they define the radial coordinates.
        !
        ! W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, Numerical Recipes in Fortran 77: The Art
        ! of Scientific Computing, 2 edition (Cambridge University Press, Cambridge England?; New York, 1992), p 246.
        !
        implicit none
        !
        integer , intent(in) :: l, m
        real(dp), intent(in) :: theta(:), phi(:)
        complex(dp), allocatable :: Ylm(:)
        integer :: n
        !
        n = size(theta)
        !
        allocate(Ylm(n))
        !
        Ylm = sqrt( real(2*l+1,dp)/fourpi * factorial(l-m)/real(factorial(l+m),dp) ) * legendre(l,m,cos(theta)) * exp(cmplx_i*m*phi)
        !
    end function  Ylm

    pure function heavi(m)
        !
        ! A. V. Podolskiy and P. Vogl, Phys. Rev. B 69, 233101 (2004). Eq 15
        !
        implicit none
        !
        integer, intent(in) :: m
        real(dp) :: heavi
        !
        if (m.ge.0) then
            heavi = 1.0_dp
        else
            heavi = 0.0_dp
        endif
        !
    end function  heavi

    ! heaviside theta functions

    pure function fermi_dirac(x) result(y)
        !
        implicit none
        !
        real(dp), intent(in) :: x
        real(dp) :: y
        real(dp) :: maxarg
        maxarg = 200.0_dp
        ! Fermi-Dirac smearing
        if (x.lt.-maxarg) then
            y = 0.0_dp
        elseif (x.gt.maxarg) then
            y = 1.0_dp
        else
            y = 1.0_dp/(1.0_dp + exp(-x))
        endif
    end function  fermi_dirac

    pure function methfessel_paxton(x,n) result(y)
        ! theta function : PRB 40, 3616 (1989).
        implicit none
        !
        real(dp), intent(in) :: x
        integer , intent(in) :: n
        real(dp) :: y

        real(dp) :: a, hp, arg, hd
        integer :: i, ni
        real(dp), parameter :: maxarg = 200.0_dp

        ! Methfessel-Paxton
        y = gauss_freq(x * sqrt(2.0_dp) )
        if (n.eq.0) return
        hd = 0.0_dp
        arg = min(maxarg, x**2)
        hp = exp(- arg)
        ni = 0
        a = 1.0_dp / sqrt (pi)
        do i=1,n
            hd = 2.0_dp*x*hp-2.0_dp*real(ni,dp)*hd
            ni = ni+1
            a = -a/(real(i,dp)*4.0_dp)
            y = y-a*hd
            hp = 2.0_dp*x*hd-2.0_dp*real(ni,dp)*hp
            ni = ni+1
        enddo
        contains
        pure function gauss_freq(x)
            !     gauss_freq(x) = (1+erf(x/sqrt(2)))/2 = erfc(-x/sqrt(2))/2
            implicit none
            !
            real(dp),intent(in) :: x
            real(dp)            :: gauss_freq
            !
            gauss_freq = 0.5_dp * erfc(-x*0.7071067811865475_dp)
            !
        end function  gauss_freq
    end function  methfessel_paxton

    pure function marzari_vanderbilt(x) result(y)
        ! theta function: PRL 82, 3296 (1999)
        ! 1/2*erf(x-1/sqrt(2)) + 1/sqrt(2*pi)*exp(-(x-1/sqrt(2))**2) + 1/2
        implicit none
        !
        real(dp), intent(in) :: x
        real(dp) :: y
        real(dp) :: arg, xp
        real(dp) :: maxarg
        maxarg = 200.0_dp
        ! Cold smearing
         xp = x - 1.0_dp / sqrt (2.0_dp)
         arg = min (maxarg, xp**2)
         y = 0.5d0 * erf (xp) + 1.0_dp / sqrt (2.0_dp * pi) * exp (- arg) + 0.5d0
    end function  marzari_vanderbilt

    pure function erf(x)
      !---------------------------------------------------------------------
      !
      !     Error function - computed from the rational approximations of
      !     W. J. Cody, Math. Comp. 22 (1969), pages 631-637.
      !
      !     for abs(x) le 0.47 erf is calculated directly
      !     for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)
      !
      implicit none
      real(dp), intent(in) :: x
      real(dp) :: x2, p1(4), q1(4)
      real(dp) :: erf
      !
      p1 = [2.426679552305318E2_dp, 2.197926161829415E1_dp, 6.996383488619136_dp,  -3.560984370181538E-2_dp]
      q1 = [2.150588758698612E2_dp, 9.116490540451490E1_dp, 1.508279763040779E1_dp, 1.000000000000000_dp]
      !
      if (abs(x).gt.6.0_dp) then
         !  erf(6)=1-10^(-17) cannot be distinguished from 1
         erf = sign (1.0_dp, x)
      else
         if (abs(x).le.0.47_dp) then
            x2 = x**2
            erf = x*(p1(1)+x2*(p1(2)+x2*(p1(3)+x2*p1(4))))/(q1(1)+x2*(q1(2)+x2*(q1(3)+x2*q1(4))))
         else
            erf = 1.0_dp - erfc(x)
         endif
      endif
    end function  erf

    pure function erfc(x)
      !     erfc(x) = 1-erf(x)  - See comments in erf
      implicit none
      real(dp),intent(in) :: x
      real(dp)            :: erfc
      real(dp) :: ax, xm2, p2(8), q2(8), p3(5), q3(5), pim1
      !
      p2 = [ 3.004592610201616E2_dp,  4.519189537118719E2_dp,  3.393208167343437E2_dp,  1.529892850469404E2_dp,  4.316222722205674E1_dp,  7.211758250883094_dp,    5.641955174789740E-1_dp,-1.368648573827167E-7_dp]
      q2 = [ 3.004592609569833E2_dp,  7.909509253278980E2_dp,  9.313540948506096E2_dp,  6.389802644656312E2_dp,  2.775854447439876E2_dp,  7.700015293522947E1_dp,  1.278272731962942E1_dp,  1.000000000000000_dp]
      p3 = [-2.996107077035422E-3_dp,-4.947309106232507E-2_dp,  -2.269565935396869E-1_dp,-2.786613086096478E-1_dp,  -2.231924597341847E-2_dp]
      q3 = [ 1.062092305284679E-2_dp, 1.913089261078298E-1_dp,  1.051675107067932_dp,    1.987332018171353_dp,     1.000000000000000_dp]
      pim1 = 0.56418958354775629_dp  ! sqrt(1/pi)
      ax = abs(x)
      if (ax > 26.0_dp) then
         !  erfc(26.0)=10^(-296); erfc(9.0)=10^(-37);
         erfc = 0.0_dp
      elseif (ax.gt.4.0_dp) then
         xm2=(1.0_dp/ax)**2
         erfc=(1.0_dp/ax)*exp(-2.0_dp)*(pim1+xm2*(p3(1)+xm2*(p3(2)+xm2*(p3(3)+xm2*(p3(4)+xm2*p3(5)))))/(q3(1)+xm2*(q3(2)+xm2*(q3(3)+xm2*(q3(4)+xm2*q3(5))))))
      elseif(ax.gt.0.47_dp)then
         erfc=exp(-x**2)*(p2(1)+ax*(p2(2)+ax*(p2(3)+ax*(p2(4)+ax*(p2(5)+ax*(p2(6)+ax*(p2(7)+ax*p2(8))))))))/(q2(1)+ax*(q2(2)+ax*(q2(3)+ax*(q2(4)+ax*(q2(5)+ax*(q2(6)+ax*(q2(7)+ax*q2(8))))))))
      else
         erfc=1.0_dp-erf(ax)
      endif
      ! erf(-x)=-erf(x)  =>  erfc(-x) = 2-erfc(x)
      if (x < 0.0_dp) erfc = 2.0_dp - erfc
    end function  erfc

    ! delta functions and theta function derivatives

    pure function fermi_dirac_dydx(x) result(y)
        ! derivative of Fermi-Dirac function: 0.5/(1.0+cosh(x))
        implicit none
        real(dp), intent(in) :: x
        real(dp) :: y
        !
        if (abs(x).le.36.0_dp) then
            y = 1.0_dp/(2.0_dp+exp(-x)+exp(+x))
        else
            y = 0.0_dp
        endif
    end function  fermi_dirac_dydx

    pure function marzari_vanderbilt_dydx(x) result(y)
        ! 1/sqrt(pi)*exp(-(x-1/sqrt(2))**2)*(2-sqrt(2)*x)
        implicit none
        real(dp), intent(in) :: x
        real(dp) :: y
        real(dp) :: arg
        real(dp) :: sqrtpm1
        !
        sqrtpm1 = 0.564189583547756_dp ! 1/sqrt(pi)
        arg = min(200.0_dp,(x-1.0_dp/sqrt(2.0_dp))**2)
        y = sqrtpm1*exp(-arg)*(2.0_dp-sqrt(2.0_dp)*x)
    end function  marzari_vanderbilt_dydx

    pure function methfessel_paxton_dydx(x,n) result(y)
        ! derivative of the corresponding Methfessel-Paxton wgauss
        implicit none
        real(dp), intent(in) :: x
        integer , intent(in) :: n
        real(dp) :: y
        real(dp) :: a, arg, hp, hd
        integer :: i, ni
        real(dp) :: sqrtpm1
        !
        sqrtpm1 = 0.564189583547756_dp ! 1/sqrt(pi)
        arg = min(200.0_dp, x**2)
        y = exp(-arg)*sqrtpm1
        if (n.eq.0) return
        hd = 0.0_dp
        hp = exp(-arg)
        ni = 0
        a  = sqrtpm1
        do i = 1, n
            hd = 2.0_dp*x*hp-2.0_dp*real(ni,dp)*hd
            ni = ni+1
            a  = -a/(real(i,dp)*4.0_dp)
            hp = 2.0_dp*x*hd-2.0_dp*real(ni,dp)*hp
            ni = ni+1
            y  = y+a*hp
        enddo
    end function  methfessel_paxton_dydx

    pure function gauss(x) result(y)
        ! derivative of Fermi-Dirac function: exp( - (x-x_o)^2 / (2*sigma^2) ) / sigma * sqrt( 2 * pi )
        ! degauss = (sqrt(2)*sigma) => exp(-(x-xo)**2/degauss**2)/(degauss*sqrt(pi))
        ! set degauss = 1 and center xo = 0
        implicit none
        real(dp), intent(in) :: x
        real(dp) :: y
        !
        if (abs(x).le.7.0_dp) then
            y = exp(-x**2.0_dp)/sqrtpi
            ! in order to avoid problems for large values of x in the e
        else
            y = 0.0_dp
        endif
    end function  gauss

    pure function lorentz(x) result(y)
        ! x = ( E(j,k)-Ep(i) )/degauss
        ! derivative of Fermi-Dirac function: 1/( ((x-xo)/degauss)**2 + 1) * 1/(pi * degauss),
        ! notice that the 1/degauss factor is taken outside of the equation
        implicit none
        real(dp), intent(in) :: x
        real(dp) :: y
        !
        y = 1.0_dp/( pi * ( x**2 + 1 ) )
        !
    end function  lorentz

    ! $$$$$$\ $$\   $$\ $$$$$$$$\ $$$$$$$$\  $$$$$$\  $$$$$$$\   $$$$$$\ $$$$$$$$\ $$$$$$\  $$$$$$\  $$\   $$\
    !  \_$$  _|$$$\  $$ |\__$$  __|$$  _____|$$  __$$\ $$  __$$\ $$  __$$\\__$$  __|\_$$  _|$$  __$$\ $$$\  $$ |
    !    $$ |  $$$$\ $$ |   $$ |   $$ |      $$ /  \__|$$ |  $$ |$$ /  $$ |  $$ |     $$ |  $$ /  $$ |$$$$\ $$ |
    !    $$ |  $$ $$\$$ |   $$ |   $$$$$\    $$ |$$$$\ $$$$$$$  |$$$$$$$$ |  $$ |     $$ |  $$ |  $$ |$$ $$\$$ |
    !    $$ |  $$ \$$$$ |   $$ |   $$  __|   $$ |\_$$ |$$  __$$< $$  __$$ |  $$ |     $$ |  $$ |  $$ |$$ \$$$$ |
    !    $$ |  $$ |\$$$ |   $$ |   $$ |      $$ |  $$ |$$ |  $$ |$$ |  $$ |  $$ |     $$ |  $$ |  $$ |$$ |\$$$ |
    !  $$$$$$\ $$ | \$$ |   $$ |   $$$$$$$$\ \$$$$$$  |$$ |  $$ |$$ |  $$ |  $$ |   $$$$$$\  $$$$$$  |$$ | \$$ |
    !  \______|\__|  \__|   \__|   \________| \______/ \__|  \__|\__|  \__|  \__|   \______| \______/ \__|  \__|

    ! gaussian method

    function       get_dos_quick(Ep,E,kptw,degauss,flags) result(D)
        ! flags = fermi, mp, mv, gauss, lorentz
        implicit none
        !
        real(dp), intent(in) :: Ep(:)       ! probing energies
        real(dp), intent(in) :: E(:,:)      ! E(nbands,nkpts) band energies
        real(dp), intent(in) :: kptw(:)     ! kptw(nkpts) normalized kpoint weights
        real(dp), intent(in) :: degauss
        character(*), intent(in) :: flags
        real(dp), allocatable :: D(:)
        integer  :: nEs
        integer  :: nbands
        integer  :: nkpts
        integer  :: i,j,k
        real(dp) :: xp
        ! get number of probing energies
        nEs = size(Ep)
        ! get n of band
        nbands = size(E,1)
        ! get n of kpoints
        nkpts = size(kptw)
        ! allocate space for dos
        allocate(D(nEs))
        !$OMP PARALLEL PRIVATE(i,j,k) SHARED(D,E,Ep,nEs,kptw)
        !$OMP DO
        do i = 1, nEs
            ! initialize D(i)
            D(i) = 0.0_dp
            do j = 1, nbands
            do k = 1, nkpts
                xp = ( E(j,k)-Ep(i) )/degauss
                ! get contribution
                if     (index(flags,'fermi').ne.0) then
                    D(i) = D(i) + kptw(k) * fermi_dirac_dydx(x=xp)
                elseif (index(flags,'mp').ne.0) then
                    D(i) = D(i) + kptw(k) * methfessel_paxton_dydx(x=xp,n=1)
                elseif (index(flags,'mv').ne.0) then
                    D(i) = D(i) + kptw(k) * marzari_vanderbilt_dydx(x=xp)
                elseif (index(flags,'gauss').ne.0) then
                    D(i) = D(i) + kptw(k) * gauss(x=xp)
                elseif (index(flags,'lorentz').ne.0) then
                    D(i) = D(i) + kptw(k) * lorentz(x=xp)
                else
                    stop 'ERROR [get_dos_quick]: flag not recognized'
                endif
            enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        D = D / degauss
    end function   get_dos_quick

    function       get_pdos_quick(Ep,E,kptw,degauss,W,flags) result(pD)
        ! flags = fermi, mp, mv, gauss, lorentz
        implicit none
        !
        real(dp), intent(in) :: Ep(:)    ! probing energies
        real(dp), intent(in) :: E(:,:)   ! E(nbands,nkpts) band energies
        real(dp), intent(in) :: W(:,:,:) ! W(nprojections,nbands,nkpts) band weights
        real(dp), intent(in) :: kptw(:)  ! kptw(nkpts) normalized kpoint weights
        real(dp), intent(in) :: degauss
        character(*), intent(in) :: flags
        real(dp), allocatable :: pD(:,:)
        integer  :: nEs
        integer  :: nbands
        integer  :: nkpts
        integer  :: nprojections
        integer  :: i,j,k,l
        real(dp) :: xp,yp
        ! get number of probing energies
        nEs = size(Ep)
        ! get n of band
        nbands = size(E,1)
        ! get n of kpoints
        nkpts = size(kptw)
        ! numbe of projections
        nprojections = size(W,1)
        ! allocate space for dos
        allocate(pD(nprojections,nEs))
        ! initialize
        pD = 0.0_dp
        !$OMP PARALLEL PRIVATE(i,j,k,l) SHARED(pD,E,Ep,nEs,kptw)
        !$OMP DO
        do i = 1, nEs
            do j = 1, nbands
            do k = 1, nkpts
                xp = ( E(j,k)-Ep(i) )/degauss
                ! get contribution from this band
                if     (index(flags,'fermi').ne.0) then
                    yp = fermi_dirac_dydx(x=xp)
                elseif (index(flags,'mp').ne.0) then
                    yp = methfessel_paxton_dydx(x=xp,n=1)
                elseif (index(flags,'mv').ne.0) then
                    yp = marzari_vanderbilt_dydx(x=xp)
                elseif (index(flags,'gauss').ne.0) then
                    yp = gauss(x=xp)
                elseif (index(flags,'lorentz').ne.0) then
                    yp = lorentz(x=xp)
                else
                    stop 'ERROR [get_dos_quick]: flag not recognized'
                endif
                ! multiply by kpoint weight
                yp = yp * kptw(k)
                ! multiply by band weights
                do l = 1, nprojections
                    pD(l,i) = pD(l,i) + yp * W(l,j,k)
                enddo
            enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        pD = pD / degauss
    end function   get_pdos_quick

    ! DOS properties

    function       get_Ef(E,tet,tetw,nelecs) result(Ef)
        ! adapted from qe
        implicit none
        !
        real(dp), intent(in) :: E(:,:)
        integer , intent(in) :: tet(:,:)
        real(dp), intent(in) :: tetw(:)
        real(dp), intent(in) :: nelecs
        real(dp) :: Ef
        integer  :: i
        integer  :: ntets
        integer  :: nbands
        integer  :: nkpts
        real(dp) :: Emin    ! absolute lowest band energy
        real(dp) :: Emax    ! absolute highest band energy
        real(dp) :: elw     ! lower bound for EF
        real(dp) :: eup     ! upper bound for EF
        real(dp) :: idos_up ! number of states for Ef=eup
        real(dp) :: idos_lw ! number of states for Ef=elw
        real(dp) :: idos_ep ! number of states for Ef=(elw+eup)/2
        real(dp) :: try_Ef  ! new Ef to try
        real(dp) :: crit    ! stopping criterion
        ! get bands, kpts, tetrahedra
        nbands = size(E,1)
        nkpts  = size(E,2)
        ntets  = sizE(tet,2)
        ! get lower and upper bracketing bounds
        Emin = minval(E)
        Emax = maxval(E)
        elw = Emin
        eup = Emax
        ! initailize bracketing
        idos_up = get_tet_idos_engine(Ep=eup,E=E,tet=tet,tetw=tetw)
        idos_lw = get_tet_idos_engine(Ep=elw,E=E,tet=tet,tetw=tetw)
        crit = 1.0d+10
        ! search for Ef by gradualling reducing bracketted region
        search : do i = 1, 1000 ! hard-coded 1000 maximum iterations
            try_Ef = (eup+elw) / 2.0_dp
            idos_ep = get_tet_idos_engine(Ep=try_Ef,E=E,tet=tet,tetw=tetw)
            if (abs(idos_ep-nelecs).lt.crit) then
                crit = abs(idos_ep-nelecs)
                Ef = try_Ef
            endif
            ! converged
            if (abs(idos_ep-nelecs).lt.tiny) then
                exit search
            elseif((idos_ep-nelecs).lt.-tiny) then
                elw = try_Ef
            else
                eup = try_Ef
            endif
        enddo search
        ! check convergence
        if (abs(idos_ep-nelecs).gt.tiny) stop 'ERROR [XX] failed to converge on Ef.'
        ! check that Ef < highest band
        if (Ef.gt.Emax) stop 'ERROR [XX]: Ef is above the highest band energy'
        !
    end function   get_Ef

    ! tetrahedra methods

    function       get_dos_tet(Ep,E,tet,tetw) result(D)
        !
        implicit none
        !
        real(dp), intent(in) :: Ep(:) ! probing energies
        real(dp), intent(in) :: E(:,:) ! E(nbands,nkpts) band energies
        integer , intent(in) :: tet(:,:) ! tet(:,ntets) tetrahedra conenctivity
        real(dp), intent(in) :: tetw(:) ! tetw(ntets) weight of each tetrahedron
        real(dp), allocatable :: D(:)
        integer :: nEs
        integer :: i
        ! get number of probing energies
        nEs = size(Ep)
        ! allocate space for dos
        allocate(D(nEs))
        !$OMP PARALLEL PRIVATE(i) SHARED(D,nEs,E,Ep,tet,tetw)
        !$OMP DO
        do i = 1, nEs
            ! get dos
            D(i) = get_tet_dos_engine(Ep=Ep(i),E=E,tet=tet,tetw=tetw)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    end function   get_dos_tet

    function       get_idos_tet(Ep,E,tet,tetw) result(iD)
        !
        implicit none
        !
        real(dp), intent(in) :: Ep(:) ! probing energies
        real(dp), intent(in) :: E(:,:) ! E(nbands,nkpts) band energies
        integer , intent(in) :: tet(:,:) ! tet(:,ntets) tetrahedra conenctivity
        real(dp), intent(in) :: tetw(:) ! tetw(ntets) weight of each tetrahedron
        real(dp), allocatable :: iD(:)
        integer :: nEs
        integer :: i
        ! get number of probing energies
        nEs = size(Ep)
        ! allocate space for dos
        allocate(iD(nEs))
        !$OMP PARALLEL PRIVATE(i) SHARED(iD,nEs,E,Ep,tet,tetw)
        !$OMP DO
        do i = 1, nEs
            ! get dos
            iD(i) = get_tet_idos_engine(Ep=Ep(i),E=E,tet=tet,tetw=tetw)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    end function   get_idos_tet

    function       get_pdos_tet(Ep,E,tet,tetw,weight) result(pD)
        !
        implicit none
        !
        real(dp), intent(in) :: Ep(:) ! probing energies
        real(dp), intent(in) :: E(:,:) ! E(nbands,nkpts) band energies
        integer , intent(in) :: tet(:,:) ! tet(:,ntets) tetrahedra conenctivity
        real(dp), intent(in) :: tetw(:) ! tetw(ntets) weight of each tetrahedron
        real(dp), intent(in) :: weight(:,:,:) ! weights(nprojections,nbands,nkpts), weights (per band per kpoint) for projected dos: can be absolute square of TB or BvK eigenvectors
        real(dp), allocatable :: pD(:,:)
        integer  :: nEs
        integer  :: nbands
        integer  :: ntets
        integer  :: nprojections
        real(dp) :: wc(4) ! corner tetrahedron weights
        integer  :: i,j,k,m
        ! get number of probing energies
        nEs = size(Ep)
        ! get n of band
        nbands = size(E,1)
        ! get number of tetrahedra
        ntets = size(tet,2)
        ! get number of projections
        nprojections = size(weight,1)
        ! allocate space
        allocate(pD(nprojections,nEs))
        ! loop over energies, bands, tetrahedra
        !$OMP PARALLEL PRIVATE(i,j,k,m,wc) SHARED(pD,E,Ep,nEs,tetw,weight)
        !$OMP DO
        do i = 1, nEs
            ! initialize D(i)
            pD(:,i) = 0.0_dp
            do j = 1, nbands
            do k = 1, ntets
                ! get tetrahedron corner weights
                wc = get_delta_wc(Ep=Ep(i),Ec=E(j,tet(:,k)))
                ! loop over projections, increment pDOS with contributions from band j in tetrahedron k
                do m = 1, nprojections
                    pD(m,i) = pD(m,i) + tetw(k) * sum(wc * weight(m,j,tet(:,k)) )
                enddo
            enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    end function   get_pdos_tet

    pure function  get_theta_wc(Ep,Ec,ntets) result(wc)
        ! get linear tetrahedron corner weights for delta function with Blochl corrections
        implicit none
        !
        real(dp), intent(in) :: Ep    ! probing energy
        real(dp), intent(in) :: Ec(4) ! corner energies
        integer , intent(in) :: ntets
        real(dp) :: wc(4)
        real(dp) :: dosEp
        real(dp) :: etot
        real(dp) :: c1,c2,c3,c4
        real(dp) :: e1,e2,e3,e4
        integer  :: k1,k2,k3,k4
        integer  :: inds(4)
        ! rank corner weights in increasing order
        inds = fourrank(Ec)
        ! sort weights
        e1 = Ec(inds(1))
        e2 = Ec(inds(2))
        e3 = Ec(inds(3))
        e4 = Ec(inds(4))
        ! k1-k4 are the irreducible k-points corresponding to e1-e4
        k1 = inds(1)
        k2 = inds(2)
        k3 = inds(3)
        k4 = inds(4)
        ! calculate weights wc
        if     (Ep.le.e1) then
            ! Eq B1 from Blochl's PhysRevB.49.16223
            wc(k1)=0.0_dp
            wc(k2)=0.0_dp
            wc(k3)=0.0_dp
            wc(k4)=0.0_dp
        elseif (Ep.le.e2) then
            ! Eq B6 from Blochl's PhysRevB.49.16223
            c4=0.25_dp/ntets*(Ep-e1)**3/(e2-e1)/(e3-e1)/(e4-e1)
            ! Eq C2 from Blochl's PhysRevB.49.16223
            dosEp=3.0_dp/ntets*(Ep-e1)**2/(e2-e1)/(e3-e1)/(e4-e1)
            ! shortcut for Blochl corrections (Eq.22 PhysRevB.49.16223)
            etot=e1+e2+e3+e4
            ! Eq B2 from Blochl's PhysRevB.49.16223
            wc(k1)=c4*(4.0_dp-(Ep-e1)*(1.0_dp/(e2-e1)+1.0_dp/(e3-e1)+1.0_dp/(e4-e1)))+dosEp*(etot-4.0_dp*e1)*0.025_dp
            ! Eq B3 from Blochl's PhysRevB.49.16223
            wc(k2)=c4*(Ep-e1)/(e2-e1)+dosEp*(etot-4.0_dp*e2)*0.025_dp
            ! Eq B4 from Blochl's PhysRevB.49.16223
            wc(k3)=c4*(Ep-e1)/(e3-e1)+dosEp*(etot-4.0_dp*e3)*0.025_dp
            ! Eq B5 from Blochl's PhysRevB.49.16223
            wc(k4)=c4*(Ep-e1)/(e4-e1)+dosEp*(etot-4.0_dp*e4)*0.025_dp
        elseif (Ep.le.e3) then
            ! Eq B11 from Blochl's PhysRevB.49.16223
            c1=0.25_dp/ntets*(Ep-e1)**2/(e4-e1)/(e3-e1)
            ! Eq B12 from Blochl's PhysRevB.49.16223
            c2=0.25_dp/ntets*(Ep-e1)*(Ep-e2)*(e3-Ep)/(e4-e1)/(e3-e2)/(e3-e1)
            ! Eq B13 from Blochl's PhysRevB.49.16223
            c3=0.25_dp/ntets*(Ep-e2)**2*(e4-Ep)/(e4-e2)/(e3-e2)/(e4-e1)
            ! Eq C3 from Blochl's PhysRevB.49.16223
            dosEp=1.0_dp/ntets/(e3-e1)/(e4-e1)*(3.0_dp*(e2-e1)+6.0_dp*(Ep-e2)-3.0_dp*(e3-e1+e4-e2)*(Ep-e2)**2/(e3-e2)/(e4-e2))
            ! shortcut for Blochl corrections (Eq.22 PhysRevB.49.16223)
            etot=e1+e2+e3+e4
            ! Eq B7 from Blochl's PhysRevB.49.16223
            wc(k1)=c1+(c1+c2)*(e3-Ep)/(e3-e1)+(c1+c2+c3)*(e4-Ep)/(e4-e1)+dosEp*(etot-4.0_dp*e1)*0.025_dp
            ! Eq B8 from Blochl's PhysRevB.49.16223
            wc(k2)=c1+c2+c3+(c2+c3)*(e3-Ep)/(e3-e2)+c3*(e4-Ep)/(e4-e2)+dosEp*(etot-4.0_dp*e2)*0.025_dp
            ! Eq B9 from Blochl's PhysRevB.49.16223
            wc(k3)=(c1+c2)*(Ep-e1)/(e3-e1)+(c2+c3)*(Ep-e2)/(e3-e2)+dosEp*(etot-4.0_dp*e3)*0.025_dp
            ! Eq B10 from Blochl's PhysRevB.49.16223
            wc(k4)=(c1+c2+c3)*(Ep-e1)/(e4-e1)+c3*(Ep-e2)/(e4-e2)+dosEp*(etot-4.0_dp*e4)*0.025_dp
        elseif (Ep.le.e4) then
            ! Eq B18 from Blochl's PhysRevB.49.16223
            c4=0.25_dp/ntets*(e4-Ep)**3/(e4-e1)/(e4-e2)/(e4-e3)
            ! Eq C4 from Blochl's PhysRevB.49.16223
            dosEp=3.0_dp/ntets*(e4-Ep)**2/(e4-e1)/(e4-e2)/(e4-e3)
            ! shortcut for Blochl corrections (Eq.22 PhysRevB.49.16223)
            etot=e1+e2+e3+e4
            ! Eq B14 from Blochl's PhysRevB.49.16223
            wc(k1)=0.25_dp/ntets-c4*(e4-Ep)/(e4-e1)+dosEp*(etot-4.0_dp*e1)*0.025_dp
            ! Eq B15 from Blochl's PhysRevB.49.16223
            wc(k2)=0.25_dp/ntets-c4*(e4-Ep)/(e4-e2)+dosEp*(etot-4.0_dp*e2)*0.025_dp
            ! Eq B16 from Blochl's PhysRevB.49.16223
            wc(k3)=0.25_dp/ntets-c4*(e4-Ep)/(e4-e3)+dosEp*(etot-4.0_dp*e3)*0.025_dp
            ! Eq B17 from Blochl's PhysRevB.49.16223
            wc(k4)=0.25_dp/ntets-c4*(4.0_dp-(e4-Ep)*(1.0_dp/(e4-e1)+1.0_dp/(e4-e2)+1.0_dp/(e4-e3)))+dosEp*(etot-4.0_dp*e4)*0.025_dp
        elseif (Ep.ge.e4) then
            ! Eq B19 from Blochl's PhysRevB.49.16223
            wc(k1)=0.25_dp/ntets
            wc(k2)=0.25_dp/ntets
            wc(k3)=0.25_dp/ntets
            wc(k4)=0.25_dp/ntets
        endif
        !
        wc = wc * ntets
        !
    end function   get_theta_wc

    pure function  get_delta_wc(Ep,Ec) result(wc)
        ! get linear tetrahedron corner weights for delta function with Blochl corrections
        implicit none
        !
        real(dp), intent(in) :: Ep    ! probing energy
        real(dp), intent(in) :: Ec(4) ! corner energies
        real(dp) :: wc(4)
        real(dp) :: dosEp
        real(dp) :: e1,e2,e3,e4
        integer  :: k1,k2,k3,k4
        real(dp) :: o13,f12,f13,f14,f21,f23,f31,f32,f34,f24,f41,f42
        integer  :: inds(4)
        ! initialize shortcut
        o13  = 1.0_dp/3.0_dp
        ! rank corner weights in increasing order
        inds = fourrank(Ec)
        ! sort weights
        e1 = Ec(inds(1))
        e2 = Ec(inds(2))
        e3 = Ec(inds(3))
        e4 = Ec(inds(4))
        ! k1-k4 are the irreducible k-points corresponding to e1-e4
        k1 = inds(1)
        k2 = inds(2)
        k3 = inds(3)
        k4 = inds(4)
        ! calculate weights wc
        if     (Ep.le.e1) then
            ! Eq B1 from Blochl's PhysRevB.49.16223
            wc(k1)=0.0_dp
            wc(k2)=0.0_dp
            wc(k3)=0.0_dp
            wc(k4)=0.0_dp
        elseif (Ep.lt.e2) then
            f12 = (Ep-e2)/(e1-e2)
            f21 = 1.0_dp - f12
            f13 = (Ep-e3)/(e1-e3)
            f31 = 1.0_dp - f13
            f14 = (Ep-e4)/(e1-e4)
            f41 = 1.0_dp - f14
            dosEp  = 3.0_dp * f21 * f31 * f41 / (Ep-e1)
            wc(k1) = o13 * (f12 + f13 + f14)
            wc(k2) = o13 * f21
            wc(k3) = o13 * f31
            wc(k4) = o13 * f41
            wc = wc * dosEp
        elseif (Ep.lt.e3) then
            f13 = (Ep-e3)/(e1-e3)
            f31 = 1.0_dp - f13
            f14 = (Ep-e4)/(e1-e4)
            f41 = 1.0_dp-f14
            f23 = (Ep-e3)/(e2-e3)
            f32 = 1.0_dp - f23
            f24 = (Ep-e4)/(e2-e4)
            f42 = 1.0_dp - f24
            dosEp  = 3.0_dp * (f23*f31 + f32*f24)
            wc(k1) = f14 * o13 + f13*f31*f23 / dosEp
            wc(k2) = f23 * o13 + f24*f24*f32 / dosEp
            wc(k3) = f32 * o13 + f31*f31*f23 / dosEp
            wc(k4) = f41 * o13 + f42*f24*f32 / dosEp
            dosEp  = dosEp / (e4-e1)
            wc = wc * dosEp
        elseif (Ep.lt.e4) then
            f14 = (Ep-e4)/(e1-e4)
            f24 = (Ep-e4)/(e2-e4)
            f34 = (Ep-e4)/(e3-e4)
            dosEp  = 3.0_dp * f14 * f24 * f34 / (e4-Ep)
            wc(k1) = f14 * o13
            wc(k2) = f24 * o13
            wc(k3) = f34 * o13
            wc(k4) = (3.0_dp - f14 - f24 - f34 ) * o13
            wc = wc * dosEp
        elseif (Ep.gt.e4) then
            wc(k1)=0.0_dp
            wc(k2)=0.0_dp
            wc(k3)=0.0_dp
            wc(k4)=0.0_dp
       endif
       !
    end function   get_delta_wc

    pure function  get_tet_idos_engine(Ep,E,tet,tetw) result(idos)
        ! return integrated DOS at Ep
        implicit none
        !
        real(dp), intent(in) :: Ep
        real(dp), intent(in) :: E(:,:)
        integer , intent(in) :: tet(:,:)
        real(dp), intent(in) :: tetw(:) ! weight of irreducible tetrahedron
        real(dp) :: idos
        real(dp) :: Ec(4)
        real(dp) :: e1,e2,e3,e4
        integer  :: i,j
        integer  :: nbands
        integer  :: ntets
        ! get bands, kpts, tetrahedra
        nbands = size(E,1)
        ntets  = sizE(tet,2)
        !
        idos = 0.0_dp
        do j = 1, ntets
            do i = 1, nbands
                ! get energies at the vertices of the j-th tethedron in ascending order
                Ec(1:4) = foursort(E(i,tet(:,j)))
                ! map to
                e1 = Ec(1)
                e2 = Ec(2)
                e3 = Ec(3)
                e4 = Ec(4)
                ! calculate sum over k of the integrated charge
                if     (Ep.le.e1) then
                    ! Eq A1 from Blochl's PhysRevB.49.16223
                    ! do nothing
                elseif (Ep.le.e2) then
                    ! Eq A2 from Blochl's PhysRevB.49.16223
                    idos = idos + tetw(j)*(Ep-e1)**3/(e2-e1)/(e3-e1)/(e4-e1)
                elseif (Ep.le.e3) then
                    ! Eq A3 from Blochl's PhysRevB.49.16223
                    idos = idos + tetw(j)/(e3-e1)/(e4-e1)*((e2-e1)**2+3.0_dp*(e2-e1)*(Ep-e2)+3.0_dp*(Ep-e2)**2.0_dp-(e3-e1+e4-e2)/(e3-e2)/(e4-e2)*(Ep-e2)**3.0_dp)
                elseif (Ep.le.e4) then
                    ! Eq A4 from Blochl's PhysRevB.49.16223
                    idos = idos + tetw(j)*(1.0_dp-(e4-Ep)**3.0_dp/(e4-e1)/(e4-e2)/(e4-e3))
                elseif (Ep.ge.e4) then
                    ! Eq A5 from Blochl's PhysRevB.49.16223
                    idos = idos + tetw(j)
                endif
            enddo
        enddo
    end function   get_tet_idos_engine

    pure function  get_tet_dos_engine(Ep,E,tet,tetw) result(dosEp)
        ! returns DOS at Ep
        implicit none
        !
        real(dp), intent(in) :: Ep
        real(dp), intent(in) :: E(:,:)
        integer , intent(in) :: tet(:,:)
        real(dp), intent(in) :: tetw(:) ! weight of irreducible tetrahedron
        real(dp) :: dosEp
        real(dp) :: Ec(4)
        real(dp) :: e1,e2,e3,e4
        integer :: i,j
        integer :: nbands
        integer :: ntets
        ! get bands, kpts, tetrahedra
        nbands = size(E,1)
        ntets  = sizE(tet,2)
        !
        dosEp = 0.0_dp
        do j = 1, ntets
            do i = 1, nbands
                ! get energies at the verticies of the j-th tethedron in ascending order
                Ec(1:4) = foursort(E(i,tet(:,j)))
                ! map to
                e1 = Ec(1)
                e2 = Ec(2)
                e3 = Ec(3)
                e4 = Ec(4)
                ! calculate sum over k of the integrated charge
                if     (Ep.le.e1) then
                    ! Eq C1 from Blochl's PhysRevB.49.16223
                    dosEp = dosEp*1.0_dp
                    ! do nothing
                elseif (Ep.le.e2) then
                    ! Eq C2 from Blochl's PhysRevB.49.16223
                    dosEp = dosEp + tetw(j)*3.0_dp*(Ep-e1)**2.0_dp/(e2-e1)/(e3-e1)/(e4-e1)
                elseif (Ep.le.e3) then
                    ! Eq C3 from Blochl's PhysRevB.49.16223
                    dosEp = dosEp + tetw(j)/(e3-e1)/(e4-e1)*(3.0_dp*(e2-e1)+6.0_dp*(Ep-e2)-3.0_dp*(e3-e1+e4-e2)/(e3-e2)/(e4-e2)*(Ep-e2)**2.0_dp)
                elseif (Ep.le.e4) then
                    ! Eq C4 from Blochl's PhysRevB.49.16223
                    dosEp = dosEp + tetw(j)*(3.0_dp*(e4-Ep)**2.0_dp/(e4-e1)/(e4-e2)/(e4-e3))
                elseif (Ep.ge.e4) then
                    ! Eq C1 from Blochl's PhysRevB.49.16223
                    dosEp = dosEp*1.0_dp
                    ! do nothing
                endif
            enddo
        enddo
        !
    end function   get_tet_dos_engine

    ! $$\      $$\ $$$$$$\  $$$$$$\   $$$$$$\  $$$$$$$$\ $$\       $$\        $$$$$$\  $$\   $$\ $$$$$$$$\  $$$$$$\  $$\   $$\  $$$$$$\
    ! $$$\    $$$ |\_$$  _|$$  __$$\ $$  __$$\ $$  _____|$$ |      $$ |      $$  __$$\ $$$\  $$ |$$  _____|$$  __$$\ $$ |  $$ |$$  __$$\
    ! $$$$\  $$$$ |  $$ |  $$ /  \__|$$ /  \__|$$ |      $$ |      $$ |      $$ /  $$ |$$$$\ $$ |$$ |      $$ /  $$ |$$ |  $$ |$$ /  \__|
    ! $$\$$\$$ $$ |  $$ |  \$$$$$$\  $$ |      $$$$$\    $$ |      $$ |      $$$$$$$$ |$$ $$\$$ |$$$$$\    $$ |  $$ |$$ |  $$ |\$$$$$$\
    ! $$ \$$$  $$ |  $$ |   \____$$\ $$ |      $$  __|   $$ |      $$ |      $$  __$$ |$$ \$$$$ |$$  __|   $$ |  $$ |$$ |  $$ | \____$$\
    ! $$ |\$  /$$ |  $$ |  $$\   $$ |$$ |  $$\ $$ |      $$ |      $$ |      $$ |  $$ |$$ |\$$$ |$$ |      $$ |  $$ |$$ |  $$ |$$\   $$ |
    ! $$ | \_/ $$ |$$$$$$\ \$$$$$$  |\$$$$$$  |$$$$$$$$\ $$$$$$$$\ $$$$$$$$\ $$ |  $$ |$$ | \$$ |$$$$$$$$\  $$$$$$  |\$$$$$$  |\$$$$$$  |
    ! \__|     \__|\______| \______/  \______/ \________|\________|\________|\__|  \__|\__|  \__|\________| \______/  \______/  \______/

    pure subroutine uc2ws(K,M,tiny)

        implicit none

        integer :: j
        real(dp), allocatable, intent(inout) :: K(:,:)
        real(dp), intent(in) :: M(3,3)
        real(dp), intent(in) :: tiny
        real(dp) :: G(3,26) , G2(26)

        ! generate mesh
        G(1,:) = real([-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1])
        G(2,:) = real([-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1])
        G(3,:) = real([-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        G = matmul(M,G)
        G2 = sum(G**2,1)

        ! perform calc
        ! do j = 1, size(K,2)
        do concurrent (j = 1:size(K,2))
            call uc2ws_engine(K(:,j),G,G2,tiny)
        enddo
    end subroutine  uc2ws

    pure subroutine uc2ws_engine(K,G,G2,tiny)

        implicit none

        real(dp), intent(inout) :: K(3)
        real(dp), intent(in) :: G(3,26), G2(26), tiny
        real(dp)  :: P
        logical :: go
        integer :: i

        go = .true.
        do while (go)
            go = .false.
            do i = 1, 26
                P = 2*(K(1)*G(1,i) + K(2)*G(2,i) + K(3)*G(3,i)) - G2(i)
                if (P .gt. 1D-8) then
                    K = K - G(:,i)
                    go = .true.
                endif
            enddo
        enddo
    end subroutine  uc2ws_engine

end module
%}

            % mex interfaces
            i=i+1; f{i} = 'uc2ws_mex'; fid=fopen([f{i},'.f90'],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid);
%{
#include "fintrf.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    use am_mex_lib , only : uc2ws

    implicit none

    mwPointer plhs(*), prhs(*)
    integer nlhs, nrhs
    mwPointer mxGetPr
    mwPointer mxCreateDoubleMatrix
    mwPointer mxGetM, mxGetN

    real*8, allocatable :: K(:,:)
    real*8 :: tiny, G(3,26) , G2(26) , M(3,3)

    ! [K] = uc2ws(K,M,tiny)
    if(nrhs .ne. 3) call mexErrMsgIdAndTxt ('MATLAB:uc2ws_mex','Three inputs required.')
    if(nlhs .gt. 1) call mexErrMsgIdAndTxt ('MATLAB:uc2ws_mex','One output produced.')
    if(mxGetM(prhs(1)) .ne. 3) call mexErrMsgIdAndTxt ('MATLAB:uc2ws_mex','Input 1 requires three-dimensional column vectors.')
    if(mxGetM(prhs(2)) .ne. 3) call mexErrMsgIdAndTxt ('MATLAB:uc2ws_mex','Input 2 requires a three-dimensional transformation matrix.')
    if(mxGetN(prhs(2)) .ne. 3) call mexErrMsgIdAndTxt ('MATLAB:uc2ws_mex','Input 2 requires a three-dimensional transformation matrix.')
    if(mxGetM(prhs(3)) .ne. 1 .or. &
     & mxGetN(prhs(3)) .ne. 1) call mexErrMsgIdAndTxt ('MATLAB:uc2ws_mex','Input 3 requires a scalar.')

    ! inputs
    ! K
    allocate( K(3,mxGetN(prhs(1))) )
    call mxCopyPtrToReal8(mxGetPr(prhs(1)),K,mxGetM(prhs(1))*mxGetN(prhs(1)))
    ! M(3,3)
    call mxCopyPtrToReal8(mxGetPr(prhs(2)),M,mxGetM(prhs(2))*mxGetN(prhs(2)))
    ! tiny
    call mxCopyPtrToReal8(mxGetPr(prhs(3)),tiny,mxGetM(prhs(3)))

    ! execute code
    call uc2ws(K,M,tiny)

    ! output results
    plhs(1) = mxCreateDoubleMatrix(mxGetM(prhs(1)),mxGetN(prhs(1)),0)
    call mxCopyReal8ToPtr(K,mxGetPr(plhs(1)),mxGetM(prhs(1))*mxGetN(prhs(1)))
end subroutine mexFunction
%}
            i=i+1; f{i} = 'get_dos_tet_mex'; fid=fopen([f{i},'.f90'],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid);
%{
#include "fintrf.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    use am_mex_lib , only : get_dos_tet

    implicit none

    mwPointer plhs(*), prhs(*)
    integer nlhs, nrhs
    mwPointer mxGetPr
    mwPointer mxCreateDoubleMatrix
    mwPointer mxGetM, mxGetN

    real*8, allocatable :: Ep(:)
    real*8, allocatable :: E(:,:)
    real*8, allocatable :: tetw(:)
    real*8, allocatable :: tet(:,:)
    real*8, allocatable :: D(:)

    integer :: i

    ! [D] = get_dos_tet(Ep,E,tet,tetw)
    if(nrhs .ne. 4) call mexErrMsgIdAndTxt ('MATLAB:get_dos_tet_mex','Four inputs required.')
    if(nlhs .gt. 1) call mexErrMsgIdAndTxt ('MATLAB:get_dos_tet_mex','One output produced.')
    if(mxGetM(prhs(3)) .ne. 4) call mexErrMsgIdAndTxt ('MATLAB:get_dos_tet_mex','Connectivity list requires four numbers in each column.')
    if(mxGetN(prhs(3)) .ne. mxGetN(prhs(4))) call mexErrMsgIdAndTxt ('MATLAB:get_dos_tet_mex','Mismatched number of tetrahedra in weights and connectivity list.')

    ! inputs
    i=0
    ! 1) Ep
    i=i+1
    allocate(   Ep( mxGetM(prhs(i))*mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),Ep, mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 2) E(nbands,nkpts)
    i=i+1
    allocate(   E( mxGetM(prhs(i)), mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),E,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 3) tet(1:4,ntets)
    i=i+1
    allocate( tet( mxGetM(prhs(i)), mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),tet,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 4) tetw(ntets)
    i=i+1
    allocate( tetw( mxGetM(prhs(i))*mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),tetw,mxGetM(prhs(i))*mxGetN(prhs(i)))

    ! execute code
    D = get_dos_tet(Ep,E,nint(tet),tetw)

    ! outputs
    i=0
    ! output results
    i=i+1
    plhs(i) = mxCreateDoubleMatrix(mxGetM(prhs(1)),mxGetN(prhs(1)),0)
    call mxCopyReal8ToPtr(D,mxGetPr(plhs(1)),mxGetM(prhs(1))*mxGetN(prhs(1)))
end subroutine mexFunction
%}
            i=i+1; f{i} = 'get_idos_tet_mex'; fid=fopen([f{i},'.f90'],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid);
%{
#include "fintrf.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    use am_mex_lib , only : get_idos_tet

    implicit none

    mwPointer plhs(*), prhs(*)
    integer nlhs, nrhs
    mwPointer mxGetPr
    mwPointer mxCreateDoubleMatrix
    mwPointer mxGetM, mxGetN

    real*8, allocatable :: Ep(:)
    real*8, allocatable :: E(:,:)
    real*8, allocatable :: tetw(:)
    real*8, allocatable :: tet(:,:)
    real*8, allocatable :: D(:)

    integer :: i

    ! [D] = get_dos_tet(Ep,E,tet,tetw)
    if(nrhs .ne. 4) call mexErrMsgIdAndTxt ('MATLAB:get_idos_tet_mex','Four inputs required.')
    if(nlhs .gt. 1) call mexErrMsgIdAndTxt ('MATLAB:get_idos_tet_mex','One output produced.')
    if(mxGetM(prhs(3)) .ne. 4) call mexErrMsgIdAndTxt ('MATLAB:get_idos_tet_mex','Connectivity list requires four numbers in each column.')
    if(mxGetN(prhs(3)) .ne. mxGetN(prhs(4))) call mexErrMsgIdAndTxt ('MATLAB:get_idos_tet_mex','Mismatched number of tetrahedra in weights and connectivity list.')

    ! inputs
    i=0
    ! 1) Ep
    i=i+1
    allocate(   Ep( mxGetM(prhs(i))*mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),Ep, mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 2) E(nbands,nkpts)
    i=i+1
    allocate(   E( mxGetM(prhs(i)), mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),E,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 3) tet(1:4,ntets)
    i=i+1
    allocate( tet( mxGetM(prhs(i)), mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),tet,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 4) tetw(ntets)
    i=i+1
    allocate( tetw( mxGetM(prhs(i))*mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),tetw,mxGetM(prhs(i))*mxGetN(prhs(i)))

    ! execute code
    D = get_idos_tet(Ep,E,nint(tet),tetw)

    ! outputs
    i=0
    ! output results
    i=i+1
    plhs(i) = mxCreateDoubleMatrix(mxGetM(prhs(1)),mxGetN(prhs(1)),0)
    call mxCopyReal8ToPtr(D,mxGetPr(plhs(1)),mxGetM(prhs(1))*mxGetN(prhs(1)))
end subroutine mexFunction
%}
            i=i+1; f{i} = 'get_pdos_tet_mex'; fid=fopen([f{i},'.f90'],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid);
%{
#include "fintrf.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    use am_mex_lib , only : get_pdos_tet

    implicit none

    mwPointer plhs(*), prhs(*)
    integer nlhs, nrhs
    mwPointer mxGetPr
    mwPointer mxCreateDoubleMatrix
    mwPointer mxGetM, mxGetN

    real*8, allocatable :: Ep(:)
    real*8, allocatable :: E(:,:)
    real*8, allocatable :: tetw(:)
    real*8, allocatable :: tet(:,:)
    real*8, allocatable :: weight(:,:,:)
    real*8, allocatable :: D(:,:)

    integer :: i

    ! [D] = get_dos_tet(Ep,E,tet,tetw)
    if(nrhs .ne. 5) call mexErrMsgIdAndTxt ('MATLAB:get_pdos_tet_mex','Five inputs required.')
    if(nlhs .gt. 1) call mexErrMsgIdAndTxt ('MATLAB:get_pdos_tet_mex','One output produced.')
    ! M x N
    if(mxGetM(prhs(3)) .ne. 4) call mexErrMsgIdAndTxt ('MATLAB:get_pdos_tet_mex','Connectivity list requires four numbers in each column.')
    if(mxGetN(prhs(3)) .ne. mxGetN(prhs(4))) call mexErrMsgIdAndTxt ('MATLAB:get_pdos_tet_mex','Mismatched number of tetrahedra in weights and connectivity list.')
    ! if(mxGetN(prhs(2)) .ne. mxGetN(prhs(5))) call mexErrMsgIdAndTxt ('MATLAB:get_pdos_tet_mex','Mismatched number of kpoints in projections weights and connectivity list.')
    
    ! inputs
    i=0
    ! 1) Ep
    i=i+1
    allocate(   Ep( mxGetM(prhs(i))*mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),Ep, mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 2) E(nbands,nkpts)
    i=i+1
    allocate(   E( mxGetM(prhs(i)), mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),E,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 3) tet(1:4,ntets)
    i=i+1
    allocate( tet( mxGetM(prhs(i)), mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),tet,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 4) tetw(ntets)
    i=i+1
    allocate( tetw( mxGetM(prhs(i))*mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),tetw,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 5) weight(nprojections,nbands,nkpts)
    i=i+1
    allocate( weight( mxGetM(prhs(i)) , mxGetM(prhs(2)) , mxGetN(prhs(2)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),weight,mxGetM(prhs(i))*mxGetM(prhs(2))*mxGetN(prhs(2)))

    ! execute code
    D = get_pdos_tet(Ep,E,nint(tet),tetw,weight)

    ! outputs
    i=0
    ! output results
    i=i+1
    plhs(i) = mxCreateDoubleMatrix(mxGetM(prhs(5)),mxGetM(prhs(1))*mxGetN(prhs(1)),0)
    call mxCopyReal8ToPtr(D,mxGetPr(plhs(1)),mxGetM(prhs(5))*mxGetM(prhs(1))*mxGetN(prhs(1)))
end subroutine mexFunction
%}

            % compile everything

            fprintf('Building mex library:\n')
                if system([FC,' ',FFLAGS,' ',DEBUG,' ',LIBS,' -c ',flib,'.f90'])
                    error(' ... %s (failed)\n',flib)
                else
                    fprintf(' ... %s (succeeded)\n',flib);
                end

            fprintf('Building mex interfaces:\n')
            for i = 1:numel(f)
                if system([FC,' ',FFLAGS,' ',DEBUG,' ',LIBS,' ',flib,'.o ',f{i},'.f90 -o ',f{i},EXT])
                    warning(' ... %s (failed)\n',f{i})
                else
                    fprintf(' ... %s (succeeded)\n',f{i});
                end
            end

        end

        function [K]   = uc2ws(K,M,flag)
            % uc2ws uses M real (reciprocal) lattice vectors to reduces K(1:3,:) vectors
            % in cartesian (reciprocal) coordinates to the definiging Wigner-Seitz cell.
            if nargin<3; flag='mex'; end
            if  contains(flag,'mex') && am_dft.usemex && exist('uc2ws_mex','file')==3
                % use mex function if available
                K = uc2ws_mex(K,M,am_lib.eps);
            else
                % generate mesh
                G = M * [-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1;
                         -1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1;
                         -1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1];
                G2=sum(G.^2,1);
                % call engine
                for j = 1:size(K,2); K(:,j) = uc2ws_engine(K(:,j),G,G2,am_lib.eps); end
            end

            function [K] = uc2ws_engine(K,G,G2,tiny)
                go = true;
                while go; go=false;
                    for i = 1:26
                        P = 2*K.'*G(:,i)-G2(i);
                        if P > + tiny
                            K = K - G(:,i); go=true;
                        end
                    end
                end
            end
        end

        function [D]   = get_dos_tet(Ep,E,tet,tetw,tetv,flag)
            % this is working!
            if nargin<6; flag='mex'; end
            if  contains(flag,'mex') && am_dft.usemex && exist('get_dos_tet_mex','file')==3
                D = get_dos_tet_mex(Ep,E,tet,tetw*tetv);
            else
                nEps = numel(Ep); D = zeros(nEps,1);
                h = waitbar(0,'Integrating...');
                for m = 1:nEps
                    waitbar(m./nEps,h,'Integrating...');
                    D(m) = get_dos_tet_engine(Ep(m),E,tet,tetw);
                end
                D = reshape(D,size(Ep));
                D = D * tetv;
                close(h);
            end

            function [dosEp] = get_dos_tet_engine(Ep,E,tet,tetw)
                nbands = size(E,1);
                ntets = size(tet,2);
                dosEp = 0;
                for j = 1:ntets
                    for i = 1:nbands
                        % get energies at the verticies of the j-th tethedron in ascending order
                        Ec(1:4) = sort(E(i,tet(:,j)));
                        % calculate sum over k of the integrated charge
                        if     Ep<=Ec(1)
                            % Eq C1 from Blochl's PhysRevB.49.16223
                            % do nothing
                        elseif Ep<=Ec(2)
                            % Eq C2 from Blochl's PhysRevB.49.16223
                            dosEp = dosEp + tetw(j)*3.0*(Ep-Ec(1)).^2/(Ec(2)-Ec(1))/(Ec(3)-Ec(1))/(Ec(4)-Ec(1));
                        elseif Ep<=Ec(3)
                            % Eq C3 from Blochl's PhysRevB.49.16223
                            dosEp = dosEp + tetw(j)/(Ec(3)-Ec(1))/(Ec(4)-Ec(1))*(3.0*(Ec(2)-Ec(1))+6.0*(Ep-Ec(2))-3.0*(Ec(3)-Ec(1)+Ec(4)-Ec(2))/(Ec(3)-Ec(2))/(Ec(4)-Ec(2))*(Ep-Ec(2)).^2);
                        elseif Ep<=Ec(4)
                            % Eq C4 from Blochl's PhysRevB.49.16223
                            dosEp = dosEp + tetw(j)*(3.0*(Ec(4)-Ep).^2/(Ec(4)-Ec(1))/(Ec(4)-Ec(2))/(Ec(4)-Ec(3)));
                        elseif Ep>=Ec(4)
                            % Eq C1 from Blochl's PhysRevB.49.16223
                            % do nothing
                        end
                    end
                end
            end
        end

        function [iD]  = get_idos_tet(Ep,E,tet,tetw,tetv,flag)
            % this is working!
            if nargin<7; flag='mex'; end
            if  contains(flag,'mex') && am_dft.usemex && exist('get_idos_tet_mex','file')==3
                iD = get_idos_tet_mex(Ep,E,tet,tetw*tetv);
            else
                nEps = numel(Ep); iD = zeros(nEps,1);
                h = waitbar(0,'Integrating...');
                for m = 1:nEps
                    waitbar(m./nEps,h,'Integrating...');
                    iD(m) = get_dos_tet_engine(Ep(m),E,tet,tetw);
                end
                iD = reshape(iD,size(Ep));
                iD = iD * tetv;
                close(h);
            end
            
            function [idos] = get_dos_tet_engine(Ep,E,tet,tetw)
                nbands = size(E,1);
                ntets = size(tet,2);
                idos = 0;
                for j = 1:ntets
                    for i = 1:nbands
                        % get energies at the verticies of the j-th tethedron in ascending order
                        Ec(1:4) = sort(E(i,tet(:,j)));
                        % calculate sum over k of the integrated charge
                        if     Ep<=Ec(1)
                            % Eq C1 from Blochl's PhysRevB.49.16223
                            % do nothing
                        elseif Ep<=Ec(2)
                            % Eq A2 from Blochl's PhysRevB.49.16223
                            idos = idos + tetw(j)*(Ep-Ec(1)).^3/(Ec(2)-Ec(1))/(Ec(3)-Ec(1))/(Ec(4)-Ec(1));
                        elseif Ep<=Ec(3)
                            % Eq A3 from Blochl's PhysRevB.49.16223
                            idos = idos + tetw(j)/(Ec(3)-Ec(1))/(Ec(4)-Ec(1))*((Ec(2)-Ec(1)).^2+3.0*(Ec(2)-Ec(1))*(Ep-Ec(2))+3.0*(Ep-Ec(2)).^2.0-(Ec(3)-Ec(1)+Ec(4)-Ec(2))/(Ec(3)-Ec(2))/(Ec(4)-Ec(2))*(Ep-Ec(2)).^3.0);
                        elseif Ep<=Ec(4)
                            % Eq A4 from Blochl's PhysRevB.49.16223
                            idos = idos + tetw(j)*(1.0-(Ec(4)-Ep).^3.0/(Ec(4)-Ec(1))/(Ec(4)-Ec(2))/(Ec(4)-Ec(3)));
                        elseif Ep>=Ec(4)
                            % Eq A5 from Blochl's PhysRevB.49.16223
                            idos = idos + tetw(j);
                        end
                    end
                end
            end
        end
        
        function [pD]  = get_pdos_tet(Ep,E,tet,tetw,tetv,projw,flag)
            %
            if nargin<7; flag='mex'; end
            if  contains(flag,'mex') && am_dft.usemex && exist('get_pdos_tet_mex','file')==3
                % projw(nprojections,nbands,nkpts)
                if size(projw,2) ~= size(E,1); error('projections and energies: nbands mismatch'); end
                if size(projw,3) ~= size(E,2); error('projections and energies: nks mismatch'); end
                pD = get_pdos_tet_mex(Ep,E,tet,tetw*tetv,projw);
            else
                % works well.
                % get number of probing energies
                nEps = numel(Ep);
                % get n of band
                nbands = size(E,1);
                % get number of tetrahedra
                ntets = size(tet,2);
                % get number of projections: projw(nprojections,nbands,nkpts)
                nprojections = size(projw,1);
                % initalize
                pD = zeros(nprojections,nEps);
                % define flat sum
                sum_ = @(x) sum(x(:));
                % loop over energies, bands, tetrahedra
                h = waitbar(0,'Integrating...');
                for i = 1:nEps
                waitbar(i./nEps,h,'Integrating...');
                for j = 1:nbands
                for k = 1:ntets
                    % get tetrahedron corner weights
                    wc = get_delta_wc(Ep(i),E(j,tet(:,k)));
                    % loop over projections, increment pDOS with contributions from band j in tetrahedron k
                    for m = 1:nprojections
                        pD(m,i) = pD(m,i) + tetw(k) * sum_(wc .* projw(m,j,tet(:,k)));
                    end
                end
                end
                end
                pD = pD*tetv/4;
                close(h);
            end
            function wc = get_delta_wc(Ep,Ec)
                % tested. works well
                % get linear tetrahedron corner weights for delta function with Blochl corrections
                % rank corner weights in increasing order
                inds = am_lib.rank4_(Ec);
                % sort weights
                e1 = Ec(inds(1)); e2 = Ec(inds(2)); e3 = Ec(inds(3)); e4 = Ec(inds(4));
                % k1-k4 are the irreducible k-points corresponding to e1-e4
                k1 = inds(1);     k2 = inds(2);     k3 = inds(3);     k4 = inds(4);
                % calculate weights wc
                if Ep<e1
                    % Eq B1 from Blochl's PhysRevB.49.16223
                    wc(k1)=0; wc(k2)=0; wc(k3)=0; wc(k4)=0;
                elseif Ep<e2
                    o13 = 1/3;
                    f12 = (Ep-e2)/(e1-e2); f21 = 1 - f12;
                    f13 = (Ep-e3)/(e1-e3); f31 = 1 - f13;
                    f14 = (Ep-e4)/(e1-e4); f41 = 1 - f14;
                    dosEp  = 3 * f21 * f31 * f41 / (Ep-e1);
                    wc(k1) = o13 * (f12 + f13 + f14);
                    wc(k2) = o13 * f21;
                    wc(k3) = o13 * f31;
                    wc(k4) = o13 * f41;
                    wc = wc * dosEp;
                elseif Ep<e3
                    o13 = 1/3;
                    f13 = (Ep-e3)/(e1-e3); f31 = 1 - f13;
                    f14 = (Ep-e4)/(e1-e4); f41 = 1 - f14;
                    f23 = (Ep-e3)/(e2-e3); f32 = 1 - f23;
                    f24 = (Ep-e4)/(e2-e4); f42 = 1 - f24;
                    dosEp  = 3 * (f23*f31 + f32*f24);
                    wc(k1) = f14 * o13 + f13*f31*f23 / dosEp;
                    wc(k2) = f23 * o13 + f24*f24*f32 / dosEp;
                    wc(k3) = f32 * o13 + f31*f31*f23 / dosEp;
                    wc(k4) = f41 * o13 + f42*f24*f32 / dosEp;
                    dosEp  = dosEp / (e4-e1);
                    wc = wc * dosEp;
                elseif Ep<e4
                    o13 = 1/3;
                    f14 = (Ep-e4)/(e1-e4);
                    f24 = (Ep-e4)/(e2-e4);
                    f34 = (Ep-e4)/(e3-e4);
                    dosEp  = 3 * f14 * f24 * f34 / (e4-Ep);
                    wc(k1) = f14 * o13;
                    wc(k2) = f24 * o13;
                    wc(k3) = f34 * o13;
                    wc(k4) = (3 - f14 - f24 - f34 ) * o13;
                    wc = wc * dosEp;
                elseif Ep>=e4
                    wc(k1)=0; wc(k2)=0; wc(k3)=0; wc(k4)=0;
                end
            end
        end

        function [ipD] = get_ipdos_tet(Ep,E,tet,tetw,tetv,projw,flag)
            % this is working!
            if nargin<7; flag='mex'; end
            if  contains(flag,'mex') && am_dft.usemex && exist('get_pdos_tet_mex','file')==3
                % use mex function if available
                error('not yet implemented');
            else
                % get number of probing energies
                nEps = numel(Ep);
                % get n of band
                nbands = size(E,1);
                % get number osf tetrahedra
                ntets = size(tet,2);
                % get number of projections: projw(nprojections,nbands,nkpts)
                nprojections = size(projw,1);
                % initalize
                ipD = zeros(nprojections,nEps);
                % define flat sum
                sum_ = @(x) sum(x(:));
                % loop over energies, bands, tetrahedra
                h = waitbar(0,'Integrating...');
                for i = 1:nEps
                waitbar(i./nEps,h,'Integrating...');
                for j = 1:nbands
                for k = 1:ntets
                    % get tetrahedron corner weights
                    wc = get_theta_wc(Ep(i),E(j,tet(:,k)),ntets);
                    % loop over projections, increment pDOS with contributions from band j in tetrahedron k
                    for m = 1:nprojections
                        ipD(m,i) = ipD(m,i) + tetw(k) * sum_(wc .* projw(m,j,tet(:,k)));
                    end
                end
                end
                end
                ipD = ipD*tetv/4;
                close(h);
            end

            function wc = get_theta_wc(Ep,Ec,ntets)
                % not yet tested.
                % get linear tetrahedron corner weights for delta function with Blochl corrections
                % rank corner weights in increasing order
                inds = am_lib.rank4_(Ec);
                % sort weights
                e1 = Ec(inds(1)); e2 = Ec(inds(2)); e3 = Ec(inds(3)); e4 = Ec(inds(4));
                % k1-k4 are the irreducible k-points corresponding to e1-e4
                k1 = inds(1);     k2 = inds(2);     k3 = inds(3);     k4 = inds(4);
                % calculate weights wc
                if     Ep<e1
                    % Eq B1 from Blochl's PhysRevB.49.16223
                    wc(k1)=0; wc(k2)=0; wc(k3)=0; wc(k4)=0;
                elseif Ep<e2
                    % Eq B6 from Blochl's PhysRevB.49.16223
                    c4=0.25/ntets*(Ep-e1)^3/(e2-e1)/(e3-e1)/(e4-e1);
                    % Eq C2 from Blochl's PhysRevB.49.16223
                    dosEp=3/ntets*(Ep-e1)^2/(e2-e1)/(e3-e1)/(e4-e1);
                    % shortcut for Blochl corrections (Eq.22 PhysRevB.49.16223)
                    etot=e1+e2+e3+e4;
                    % Eq B2 from Blochl's PhysRevB.49.16223
                    wc(k1)=c4*(4-(Ep-e1)*(1/(e2-e1)+1/(e3-e1)+1/(e4-e1)))+dosEp*(etot-4*e1)*025;
                    % Eq B3 from Blochl's PhysRevB.49.16223
                    wc(k2)=c4*(Ep-e1)/(e2-e1)+dosEp*(etot-4*e2)*025;
                    % Eq B4 from Blochl's PhysRevB.49.16223
                    wc(k3)=c4*(Ep-e1)/(e3-e1)+dosEp*(etot-4*e3)*025;
                    % Eq B5 from Blochl's PhysRevB.49.16223
                    wc(k4)=c4*(Ep-e1)/(e4-e1)+dosEp*(etot-4*e4)*025;
                elseif Ep<e3
                    % Eq B11 from Blochl's PhysRevB.49.16223
                    c1=0.25/ntets*(Ep-e1)^2/(e4-e1)/(e3-e1);
                    % Eq B12 from Blochl's PhysRevB.49.16223
                    c2=0.25/ntets*(Ep-e1)*(Ep-e2)*(e3-Ep)/(e4-e1)/(e3-e2)/(e3-e1);
                    % Eq B13 from Blochl's PhysRevB.49.16223
                    c3=0.25/ntets*(Ep-e2)^2*(e4-Ep)/(e4-e2)/(e3-e2)/(e4-e1);
                    % Eq C3 from Blochl's PhysRevB.49.16223
                    dosEp=1/ntets/(e3-e1)/(e4-e1)*(3*(e2-e1)+6*(Ep-e2)-3*(e3-e1+e4-e2)*(Ep-e2)^2/(e3-e2)/(e4-e2));
                    % shortcut for Blochl corrections (Eq.22 PhysRevB.49.16223)
                    etot=e1+e2+e3+e4;
                    % Eq B7 from Blochl's PhysRevB.49.16223
                    wc(k1)=c1+(c1+c2)*(e3-Ep)/(e3-e1)+(c1+c2+c3)*(e4-Ep)/(e4-e1)+dosEp*(etot-4*e1)*025;
                    % Eq B8 from Blochl's PhysRevB.49.16223
                    wc(k2)=c1+c2+c3+(c2+c3)*(e3-Ep)/(e3-e2)+c3*(e4-Ep)/(e4-e2)+dosEp*(etot-4*e2)*025;
                    % Eq B9 from Blochl's PhysRevB.49.16223
                    wc(k3)=(c1+c2)*(Ep-e1)/(e3-e1)+(c2+c3)*(Ep-e2)/(e3-e2)+dosEp*(etot-4*e3)*025;
                    % Eq B10 from Blochl's PhysRevB.49.16223
                    wc(k4)=(c1+c2+c3)*(Ep-e1)/(e4-e1)+c3*(Ep-e2)/(e4-e2)+dosEp*(etot-4*e4)*025;
                elseif Ep<e4
                    % Eq B18 from Blochl's PhysRevB.49.16223
                    c4=0.25/ntets*(e4-Ep)^3/(e4-e1)/(e4-e2)/(e4-e3);
                    % Eq C4 from Blochl's PhysRevB.49.16223
                    dosEp=3/ntets*(e4-Ep)^2/(e4-e1)/(e4-e2)/(e4-e3);
                    % shortcut for Blochl corrections (Eq.22 PhysRevB.49.16223)
                    etot=e1+e2+e3+e4;
                    % Eq B14 from Blochl's PhysRevB.49.16223
                    wc(k1)=0.25/ntets-c4*(e4-Ep)/(e4-e1)+dosEp*(etot-4*e1)*025;
                    % Eq B15 from Blochl's PhysRevB.49.16223
                    wc(k2)=0.25/ntets-c4*(e4-Ep)/(e4-e2)+dosEp*(etot-4*e2)*025;
                    % Eq B16 from Blochl's PhysRevB.49.16223
                    wc(k3)=0.25/ntets-c4*(e4-Ep)/(e4-e3)+dosEp*(etot-4*e3)*025;
                    % Eq B17 from Blochl's PhysRevB.49.16223
                    wc(k4)=0.25/ntets-c4*(4-(e4-Ep)*(1/(e4-e1)+1/(e4-e2)+1/(e4-e3)))+dosEp*(etot-4*e4)*025;
                elseif Ep>=e4
                    % Eq B19 from Blochl's PhysRevB.49.16223
                    wc(k1)=0.25/ntets; wc(k2)=0.25/ntets; wc(k3)=0.25/ntets; wc(k4)=0.25/ntets;
                end
                %
                wc = wc * ntets;
                %
            end
        end
        
    end

end
