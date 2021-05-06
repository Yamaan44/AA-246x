defvec = zeros(4,1)
sigmavec = zeros(4,1);
for i = 1:1
    [sigmavec(i), defvec(i)] = calc_beam(0.4*1.5, 2.3, 1.5,i)
end
figure()
plot(defvec,'ro','MarkerFaceColor', 'b')
ylabel('Deflection/Span')
figure()
plot(sigmavec,'ro','MarkerFaceColor', 'b')
ylabel('Max Stress (Pa)')

%%

function [mass_empty] = empty_weight(S_ref, b,i)
    global cst
    [c, s_htail, c_htail, s_vtail, c_vtail, l_t] = size_plane(S_ref, b);
    
    %%%% Material Properties
    rho_cf = 2000; %kg/m^3, density of carbon fiber 
    rho_g_f = 0.047468; %kg/m^2 - 1.4 oz/yd^2 fiberglass - fuselage
    rho_g_t = 0.023734; %kg/m^2 - .7 oz/yd^2 fiberglass - tails
    fiber_resin_ratio = 1; %assume fiber and resin weight equal
    
    rho_f = 25.23; %oz/cu-ft, Foam density, potentially 40 
    rho_f = rho_f* 1.00115; %oz/cu-ft to kg/m^3
    
    %%%% CALCULATE FUSELAGE WEIGHT - Currently assumes a cylinder fuselage
    %assume fuselage wall thickness is 1 cm of foam
    %To account for fuselage structure, assume fuselage covered in 1.4 oz
    %fiberglass
    
    l_fus = 0.75.*b; %ballpark for fuselage length. RC airplanes typically 75% of span
    r_fus = l_fus./8./2; %assume fuselage fineness ratio of 8
    wall_thick = 0.005; %assume fuselage wall thickness is .5 cm of foam
    vol_fuse = ((pi.*r_fus^2)-(pi.*(r_fus-wall_thick)^2))*l_fus; 
    sa_fuse = (2*pi*r_fus)*l_fus;
    
    m_fuse_glass = sa_fuse* rho_g_f;
    m_fuse_glass = m_fuse_glass + m_fuse_glass/fiber_resin_ratio; %account for expoxy weight
    m_fuse =  vol_fuse*rho_f + m_fuse_glass;
    
    %%%% WING CALCS, CURRENTLY ASSUMING RECTANGULAR SECTIONS. IF TAPERED
    %%%% DESIRED, JUST ADD MORE INPUTS TO FUNCTION AND EDIT
    %%%% cr/ct
    
    c = S_ref/b;
    cr = c;
    ct = c;
    
    %Calculate Wing Spar Volume
    t_c = .08; %assume a thickness to chord for the wing
    r_o = cst.spar_ratio*t_c*c/2; %m - estimate a thickness for the spar. 
    r_i = r_o - 0.0015875; %m - assume a constant wall thickness of 1/16 inch - from Dragonplate's website
    if r_i < 0
        r_i = 0; %make sure r_i isn't negative
    end
    vol_w_spar = ((pi*r_o^2)-(pi*r_i^2))*b;
    m_w_spar = rho_cf*vol_w_spar;
    
    xvec = linspace(0,b/2,10000); %m, x positions along wing starting from root to b/2
    cvec = ((b/2 - xvec)*cr + xvec*ct)*2/b; %m, chord at each x position along wing starting from root to b/2

    %%% USE AIRFOIL READER TO CALCULATE WING VOLUME, THEN MASS
    
    %plot(xin,yin);
    a = polyarea(cst.xin,cst.yin); %m^2/m (chord)
    vol = 2*trapz(xvec,a*cvec); %m^3, total wing volume
    m_wing = vol * rho_f +  m_w_spar; %kg, wing mass (assumes solid x section)
    
    
    %%%% TAIL CALCS, CURRENTLY ASSUMING RECTANGULAR SECTIONS. IF TAPERED
    %%%% DESIRED, JUST ADD MORE INPUTS TO FUNCTION AND EDIT
    %%%% crt/ctt/crtv/cttv
    
    bt = s_htail/c_htail; %horizontal tail span
    crt = c_htail; %horizontal root chord
    ctt = c_htail; %horizontal tip chord
    
    btv = s_vtail/c_vtail; %vert tail span
    crtv = c_vtail; %vert root chord
    cttv = c_vtail; %vert tip chord
    
    %TAIL WEIGHT, HORIZONTAL
    xtvec = linspace(0,bt/2,10000); %m, x positions along tail starting from root to bt/2
    ctvec = ((bt/2 - xtvec)*crt + xtvec*ctt)*2/bt; %m, chord at each x position along tail starting from root to b/2

    %plot(xtin,ytin);
    at = polyarea(cst.xtin,cst.ytin); %m^2/m (chord)
    volt = 2*trapz(xtvec,at*ctvec); %m^3, total tail volume
    
    %Calculate horizontal tail structure - assume covered in .7 oz
    %fiberglass
    sa_ht = 2*s_htail;
    m_ht_glass = sa_ht*rho_g_t;
    m_ht_glass = m_ht_glass + m_ht_glass/fiber_resin_ratio;

    %TAIL WEIGHT, VERTICAL
    xtvvec = linspace(0,btv,10000); %m, x positions along tail starting from root to btv
    ctvvec = ((btv - xtvvec)*crtv + xtvvec*cttv)/btv; %m, chord at each x position along tail starting from root to b/2
    
    atv = polyarea(cst.xtin,cst.ytin); %m^2/m (chord)
    voltv = trapz(xtvvec,atv*ctvvec); %m^3, total tail volume
    
    %Calculate vertical tail structure - assume covered in .7 oz
    %fiberglass
    sa_vt = 2*s_vtail;
    m_vt_glass = sa_vt*rho_g_t;
    m_vt_glass = m_vt_glass + m_vt_glass/fiber_resin_ratio;

    m_htail = volt * rho_f + m_ht_glass; %kg, tail horizontal mass (assumes solid x section)
    m_vtail = voltv * rho_f + m_vt_glass; %kg, tail vertical mass (assumes solid x section)
    
    %mass of electric components
    m_elec = .200;  
    
%     m_wing_foam = vol * rho_f
%     m_fuse
%     m_wing
%     m_htail
%     m_vtail
    
    mass_empty = m_fuse + m_wing + m_htail + m_vtail + m_elec;

end

function [sigma_max, deflection_span] = calc_beam(S_ref, weight, b,i)
    global cst
    
    %Spar Dimensions
    t_c = .08; %assume a thickness to chord for the wing
    c = S_ref/b; %wing chord - m
    L_S = weight/S_ref;
    
    %circ spar geom
    r_o = cst.spar_ratio*t_c*c/2; %m - estimate a thickness for the spar. 
    r_i = r_o - 0.0015875; %m - assume a constant wall thickness of 1/16 inch - from Dragonplate's website
    
    if r_i < 0
        r_i = 0; %make sure r_i isn't negative
    end
    
    %Spar properties
    E = 70*10^9; %MPa taken from online carbon fiber spar
    Icirc1 = pi*0.25* (0.75^4 - 0.625^4)*0.0254^4; %m^4, derived from spar geometry https://dragonplate.com/carbon-fiber-roll-wrapped-twill-tube-0625-id-x-24-gloss-finish
    Icirc2 = pi*0.25* (0.625^4 - 0.5^4)*0.0254^4; %m^4, derived from spar geometry https://dragonplate.com/carbon-fiber-roll-wrapped-twill-tube-05-id-x-24-gloss-finish
    Icirc3 = pi*0.25* (0.685^4 - 0.625^4)*0.0254^4; %m^4,  https://dragonplate.com/carbon-fiber-roll-wrapped-twill-tube-0625-id-x-96-thin-wall-gloss-finish
    Isquare = 1/12 * (0.875^4 - 0.75^4)*0.0254^4; % https://dragonplate.com/carbon-fiber-roll-wrapped-twill-square-tube-075-id-x-075-id-x-24-3
    Jcirc = pi*0.25* (r_o^4 - r_i^4); %m^4, derived from spar geometry
    
    Ivec = [Icirc1 Icirc2 Icirc3 Isquare];
    I = Ivec(i);
    
    ymax = r_o; %maximum spar distance from neutral axis

    %xlim([0 b/2]);
    %ylim([0 cr]);
    
    g_load = 1; %how much load do we expect?
    xvec = linspace(0,b/2,10000); %m, x positions along wing starting from root to b/2
    w0 = (4*weight*g_load*cst.g)/(pi*b);
    wvec = w0 * sqrt(1 - (xvec/(0.5*b)).^2);%N/m, elliptical lift distribution
    %wvec = c * L_S * ones(10000,1); %N/m, rectangular lift distribution, do not use
    wvec_flip = flip(wvec);
    
    Vvec = cumtrapz(xvec,wvec_flip); %N
    Mvec = cumtrapz(xvec,Vvec); %Nm
    
    Vvec = flip(Vvec);
    Mvec = flip(Mvec);
    
    u_primevec = cumtrapz(xvec,Mvec);
    uvec = cumtrapz(xvec,u_primevec)/(E*I); %m
    umax = max(uvec);
    
    deflection_span = umax/b; %deflection to span ratio
    
    V_max = max(Vvec);
    M_max = max(Mvec);
    %T_spar = L * d; neglecting for now in calculation as negligible

    sigma_max = M_max * ymax/I; %* 10^-6; %MPa, max compressive/tensile at root
    %tau_max =T_spar * r/J * 10^-6; %MPa, max shear due to torsion

    
    plots = 1; %toggle for debugging, need root and tip chords
    
    if plots == 1
        %chord as a fcn of span
        %figure(1)
        %fplot(@(x) ((b/2 - x)*cr + x*ct)*2/b, [0 b/2])
        %figure(2)
        %fplot(@(x) ((b/2 - x)*cr + x*ct)*2/b * L_S, [0 b/2])

        figure(3)
        subplot(3,1,1)
        plot(xvec,Vvec);
        xlabel('x (m)')
        ylabel('V(x) (N)')
        subplot(3,1,2)
        plot(xvec,Mvec);
        xlabel('x (m)')
        ylabel('M(x) (Nm)')
        subplot(3,1,3)
        plot(xvec,uvec);
        xlabel('x (m)')
        ylabel('u(x) (m)')
        
        figure(4)
        plot(xvec,wvec);
        xlabel('x (m)')
        ylabel('w(x) (N/m)')
    end
end

function [n,xin,yin] = openfile(file)
    % Read airfoil profile and place in xin, yin vectors
    fid = fopen(file);
    n = fscanf(fid,'%10f',1); %length of airfoil vectors
    data = fscanf(fid,'%10f %10f \n', [2,n]);
    fclose(fid);
    xin = data(1,:);
    big=max(xin);
    xin=xin./big;
    yin = data(2,:);
    yin=yin./big;
end