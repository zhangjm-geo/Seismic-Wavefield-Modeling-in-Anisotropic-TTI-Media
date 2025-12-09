% Elastic Velocity-Stress Finite-difference Modeling with output PML in MATLAB
% Converted from C++ code provided by Zhang Jianming

% Constants
PI = 3.1415926;
dx = 10;
dz = 10;
pml = 20;
NX = 200+2*pml;
NZ = 200+2*pml;
NT = 801;
dt = 0.0005;
N = 6;
fi = (PI / 4.0);

% Allocate memory for matrices
Txx = zeros(NX, NZ);
Txx_x = zeros(NX, NZ);
Txx_z = zeros(NX, NZ);
Tzz = zeros(NX, NZ);
Tzz_x = zeros(NX, NZ);
Tzz_z = zeros(NX, NZ);
Txz = zeros(NX, NZ);
Txz_x = zeros(NX, NZ);
Txz_z = zeros(NX, NZ);
Vx = zeros(NX, NZ);
Vx_x = zeros(NX, NZ);
Vx_z = zeros(NX, NZ);
Vz = zeros(NX, NZ);
Vz_x = zeros(NX, NZ);
Vz_z = zeros(NX, NZ);
L = zeros(NX, NZ);
M = zeros(NX, NZ);
e = zeros(NX, NZ);
C11 = zeros(NX, NZ);
C33 = zeros(NX, NZ);
C44 = zeros(NX, NZ);
C66 = zeros(NX, NZ);
C13 = zeros(NX, NZ);
Rou = zeros(NX, NZ);
C_11 = zeros(NX, NZ);
C_13 = zeros(NX, NZ);
C_15 = zeros(NX, NZ);
C_33 = zeros(NX, NZ);
C_35 = zeros(NX, NZ);
C_55 = zeros(NX, NZ);
O = zeros(NX, NZ);
Vp = zeros(NX, NZ);
Eps = zeros(NX, NZ);
Del = zeros(NX, NZ);
Vs = zeros(NX, NZ);
Gam = zeros(NX, NZ);
F = zeros(NX, NZ);
data_vx = zeros(NX, NT);
data_vz = zeros(NX, NT);

% Initialize parameters
f0 = 15;
t0 = 1.2 / f0;
SX = floor(NX / 2);
SZ = floor(NZ / 2);
R = 0.001;
Vmax = 7500;
plx = pml * dx;
plz = pml * dz;
ddx = zeros(NX, NZ);
ddz = zeros(NX, NZ);

% Initialize model parameters
for i = 1:NX
    for j = 1:NZ
        Vp(i, j) = 3000;
        Eps(i, j) = 0.4;
        Del(i, j) = 0.2;
        Vs(i, j) = 2000;
        Gam(i, j) = 0.1;
        %
        %         if j > 7.0 / 10.0 * NZ
        %             Vs(i, j) = 2800;
        %         end
        %
        %         if j > 8.0 / 10.0 * NZ
        %             Vp(i, j) = 4000;
        %         end
        %
        Rou(i, j) = 3.44 * 10^3;
        C11(i, j) = Rou(i, j) * (1.0 + 2 * Eps(i, j)) * Vp(i, j)^2;
        C33(i, j) = Rou(i, j) * Vp(i, j)^2;
        C44(i, j) = Rou(i, j) * Vs(i, j)^2;
        C66(i, j) = Rou(i, j) * (1.0 + 2 * Gam(i, j)) * Vs(i, j)^2;
        F(i, j) = 1.0 - Vs(i, j)^2 / Vp(i, j)^2;
        C13(i, j) = Rou(i, j) * Vp(i, j)^2 * sqrt(F(i, j) * (F(i, j) + 2 * Del(i, j))) - Rou(i, j) * Vs(i, j)^2;

        O(i, j) = fi;
        C_11(i, j) = C11(i, j) * cos(O(i, j))^4 + C33(i, j) * sin(O(i, j))^4 + (2 * C13(i, j) + 4 * C44(i, j)) * sin(O(i, j))^2 * cos(O(i, j))^2;
        C_13(i, j) = (C11(i, j) + C33(i, j) - 4 * C44(i, j)) * sin(O(i, j))^2 * cos(O(i, j))^2 + C13(i, j) * (sin(O(i, j))^4 + cos(O(i, j))^4);
        C_15(i, j) = (C13(i, j) - C11(i, j) + 2 * C44(i, j)) * sin(O(i, j)) * cos(O(i, j))^3 - (C13(i, j) - C33(i, j) + 2 * C44(i, j)) * sin(O(i, j))^3 * cos(O(i, j));
        C_33(i, j) = C11(i, j) * sin(O(i, j))^4 + C33(i, j) * cos(O(i, j))^4 + (2 * C13(i, j) + 4 * C44(i, j)) * sin(O(i, j))^2 * cos(O(i, j))^2;
        C_35(i, j) = (C13(i, j) - C11(i, j) + 2 * C44(i, j)) * sin(O(i, j))^3 * cos(O(i, j)) - (C13(i, j) - C33(i, j) + 2 * C44(i, j)) * sin(O(i, j)) * cos(O(i, j))^3;
        C_55(i, j) = (C11(i, j) + C33(i, j) - 2 * C13(i, j)) * sin(O(i, j))^2 * cos(O(i, j))^2 + C44(i, j) * (cos(O(i, j))^2 - sin(O(i, j))^2)^2;
    end
end

% PML setup
for i = 1:NX
    for j = 1:NZ
        if i <= pml && j <= pml
            x = pml - i;
            z = pml - j;
            ddx(i, j) = -log(R) * 3 * Vmax * x^2 / (2 * plx^2);
            ddz(i, j) = -log(R) * 3 * Vmax * z^2 / (2 * plz^2);
        elseif i <= pml && j > NZ - pml
            x = pml - i;
            z = j - (NZ - pml);
            ddx(i, j) = -log(R) * 3 * Vmax * x^2 / (2 * plx^2);
            ddz(i, j) = -log(R) * 3 * Vmax * z^2 / (2 * plz^2);
        elseif i > NX - pml && j <= pml
            x = i - (NX - pml);
            z = pml - j;
            ddx(i, j) = -log(R) * 3 * Vmax * x^2 / (2 * plx^2);
            ddz(i, j) = -log(R) * 3 * Vmax * z^2 / (2 * plz^2);
        elseif i > NX - pml && j > NZ - pml
            x = i - (NX - pml);
            z = j - (NZ - pml);
            ddx(i, j) = -log(R) * 3 * Vmax * x^2 / (2 * plx^2);
            ddz(i, j) = -log(R) * 3 * Vmax * z^2 / (2 * plz^2);
        elseif i > pml && i <= NX - pml && j <= pml
            z = pml - j;
            ddx(i, j) = 0;
            ddz(i, j) = -log(R) * 3 * Vmax * z^2 / (2 * plz^2);
        elseif i > pml && i <= NX - pml && j > NZ - pml
            z = j - (NZ - pml);
            ddx(i, j) = 0;
            ddz(i, j) = -log(R) * 3 * Vmax * z^2 / (2 * plz^2);
        elseif i <= pml && j > pml && j <= NZ - pml
            x = pml - i;
            ddx(i, j) = -log(R) * 3 * Vmax * x^2 / (2 * plx^2);
            ddz(i, j) = 0;
        elseif i > NX - pml && j > pml && j <= NZ - pml
            x = i - (NX - pml);
            ddx(i, j) = -log(R) * 3 * Vmax * x^2 / (2 * plx^2);
            ddz(i, j) = 0;
        else
            ddx(i, j) = 0;
            ddz(i, j) = 0;
        end
    end
end

% Coefficients for finite difference
cof = zeros(1, N);
if N == 1
    cof(1) = 1.0;
elseif N == 2
    cof = [1.125, -0.041666667];
elseif N == 3
    cof = [1.171875, -0.065104167, 0.0046875];
elseif N == 4
    cof = [1.196289, -0.0797526, 0.009570313, -0.0006975447];
elseif N == 5
    cof = [1.2112427, -0.08972168, 0.013842773, -0.0017656599, 0.00011867947];
elseif N == 6
    cof = [1.2213364, -0.096931458, 0.017447662, -0.0029672895, 0.00035900540, -0.000021847812];
elseif N == 7
    cof = [1.2286062, -0.10238385, 0.02047677, -0.0041789327, 0.00068945355, -0.000076922503, 0.0000042365148];
elseif N == 8
    cof = [1.2340911, -0.10664985, 0.023036367, -0.0053423856, 0.0010772712, -0.00016641888, 0.000017021711, -0.00000085234642];
end

% Main time loop
for k = 1:NT
    % Update stresses
    for i = N+1:NX-N
        for j = N+1:NZ-N
            pxVx = 0; pzVx = 0; pxVz = 0; pzVz = 0;
            for in = 0:N-1
                pxVx = pxVx + cof(in+1) * (Vx(i+in+1,j) - Vx(i-in,j));
                pxVz = pxVz + cof(in+1) * (Vz(i+in,j) - Vz(i-in-1,j));
                pzVx = pzVx + cof(in+1) * (Vx(i,j+in+1) - Vx(i,j-in));
                pzVz = pzVz + cof(in+1) * (Vz(i,j+in) - Vz(i,j-in-1));
            end
            Txx_x(i,j) = 1.0 / (1 + 0.5 * dt * ddx(i,j)) * ((1 - 0.5 * ddx(i,j) * dt) * Txx_x(i,j) + C_11(i,j) * dt / dx * pxVx + C_15(i,j) * dt / dx * pxVz);
            Txx_z(i,j) = 1.0 / (1 + 0.5 * dt * ddz(i,j)) * ((1 - 0.5 * ddz(i,j) * dt) * Txx_z(i,j) + C_13(i,j) * dt / dz * pzVz + C_15(i,j) * dt / dz * pzVx);
            Tzz_x(i,j) = 1.0 / (1 + 0.5 * dt * ddx(i,j)) * ((1 - 0.5 * ddx(i,j) * dt) * Tzz_x(i,j) + C_13(i,j) * dt / dx * pxVx + C_35(i,j) * dt / dx * pxVz);
            Tzz_z(i,j) = 1.0 / (1 + 0.5 * dt * ddz(i,j)) * ((1 - 0.5 * ddz(i,j) * dt) * Tzz_z(i,j) + C_33(i,j) * dt / dz * pzVz + C_35(i,j) * dt / dz * pzVx);
            Txx(i,j) = Txx_x(i,j) + Txx_z(i,j);
            Tzz(i,j) = Tzz_x(i,j) + Tzz_z(i,j);

            % Source term
            if i == SX && j == SZ
                tmp = (PI * f0 * (k * dt - t0))^2;
                Txx_x(i,j) = Txx_x(i,j) + exp(-tmp) * (1 - 2 * tmp);
                Txx_z(i,j) = Txx_z(i,j) + exp(-tmp) * (1 - 2 * tmp);
                Tzz_x(i,j) = Txx_x(i,j);
                Tzz_z(i,j) = Txx_z(i,j);
            end

            pxVx = 0; pzVx = 0; pxVz = 0; pzVz = 0;
            for in = 0:N-1
                pxVx = pxVx + cof(in+1) * (Vx(i+in+1,j) - Vx(i-in,j));
                pxVz = pxVz + cof(in+1) * (Vz(i+in,j) - Vz(i-in-1,j));
                pzVx = pzVx + cof(in+1) * (Vx(i,j+in+1) - Vx(i,j-in));
                pzVz = pzVz + cof(in+1) * (Vz(i,j+in) - Vz(i,j-in-1));
            end
            Txz_x(i,j) = 1.0 / (1 + 0.5 * dt * ddx(i,j)) * ((1 - 0.5 * ddx(i,j) * dt) * Txz_x(i,j) + C_15(i,j) * dt / dx * pxVx + C_55(i,j) * dt / dx * pxVz);
            Txz_z(i,j) = 1.0 / (1 + 0.5 * dt * ddz(i,j)) * ((1 - 0.5 * ddz(i,j) * dt) * Txz_z(i,j) + C_35(i,j) * dt / dz * pzVz + C_55(i,j) * dt / dz * pzVx);
            Txz(i,j) = Txz_x(i,j) + Txz_z(i,j);
        end
    end

    % Update velocities
    for i = N+1:NX-N
        for j = N+1:NZ-N
            pxTxx = 0; pzTxz = 0;
            for in = 0:N-1
                pxTxx = pxTxx + cof(in+1) * (Txx(i+in,j) - Txx(i-in-1,j));
                pzTxz = pzTxz + cof(in+1) * (Txz(i,j+in) - Txz(i,j-in-1));
            end
            Vx_x(i,j) = 1.0 / (1 + 0.5 * dt * ddx(i,j)) * ((1 - 0.5 * ddx(i,j) * dt) * Vx_x(i,j) + dt / Rou(i,j) / dx * pxTxx);
            Vx_z(i,j) = 1.0 / (1 + 0.5 * dt * ddz(i,j)) * ((1 - 0.5 * ddz(i,j) * dt) * Vx_z(i,j) + dt / Rou(i,j) / dz * pzTxz);
            Vx(i,j) = Vx_x(i,j) + Vx_z(i,j);

            pxTxz = 0; pzTzz = 0;
            for in = 0:N-1
                pxTxz = pxTxz + cof(in+1) * (Txz(i+in+1,j) - Txz(i-in,j));
                pzTzz = pzTzz + cof(in+1) * (Tzz(i,j+in+1) - Tzz(i,j-in));
            end
            Vz_x(i,j) = 1.0 / (1 + 0.5 * dt * ddx(i,j)) * ((1 - 0.5 * ddx(i,j) * dt) * Vz_x(i,j) + dt / Rou(i,j) / dx * pxTxz);
            Vz_z(i,j) = 1.0 / (1 + 0.5 * dt * ddz(i,j)) * ((1 - 0.5 * ddz(i,j) * dt) * Vz_z(i,j) + dt / Rou(i,j) / dz * pzTzz);
            Vz(i,j) = Vz_x(i,j) + Vz_z(i,j);
        end
    end



    % Plot the wavefield for animation
    ndt = 10;
    if  k>=200 & mod(k, ndt) == 0
        TT=Txx(pml+1:NX-pml,pml+1:NZ-pml)+Tzz(pml+1:NX-pml,pml+1:NZ-pml);
        imagesc(TT);
        colormap('jet');
        colorbar;
        %     caxis([-1, 1]);
        title(['Wavefield at time step: ', num2str(k)]);
        xlabel('X');
        ylabel('Z');
        drawnow;


        frame = getframe(gcf);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if k == 200
            imwrite(imind, cm, 'Wavefield.gif', 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
        else
            imwrite(imind, cm, 'Wavefield.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
    end


    
    % Record data
    for i = 1:NX
        rz = SZ;
        data_vx(i,k) = Vx(i,rz);
        data_vz(i,k) = Vz(i,rz);
    end

%     % Save data every 1000 time steps
%     if mod(k, 1000) == 0
%         Vx_filename = sprintf('Vx_%d.dat', k);
%         Vz_filename = sprintf('Vz_%d.dat', k);
%         save(Vx_filename, 'Vx', '-ascii');
%         save(Vz_filename, 'Vz', '-ascii');
%     end
end


% % Save final data
% save('Txx.dat', 'Txx', '-ascii');
% save('Tzz.dat', 'Tzz', '-ascii');
% T = 0.5 * (Txx + Tzz);
% save('T.dat', 'T', '-ascii');
% save('Txz.dat', 'Txz', '-ascii');
% save('Vx.dat', 'Vx', '-ascii');
% save('Vz.dat', 'Vz', '-ascii');
% save('data_vx.dat', 'data_vx', '-ascii');
% save('data_vz.dat', 'data_vz', '-ascii');
% save('Vp.dat', 'Vp', '-ascii');

% % Plot the wavefield
% imagesc(Vx);
% colormap('jet');
% colorbar;
% title('Vx Wavefield');
% xlabel('X');
% ylabel('Z');

