% ************************** BEGIN DESCRIPTION ****************************
% MATLAB script that generates coordinates for a coarse-grained molecular 
% dynamics model of multilayered graphene oxide in the LAMMPS "full" atom 
% style (last updated 7/24/25 by Heather White)
%
% Supervised by Sinan Keten (s-keten@u.northwestern.edu)
% Developed by Luis Ruiz (lar181@miami.edu), Zhaoxu Meng
% (zmeng@clemson.edu, and Heather White (heatherosa37@gmail.com)
%
% Select publications (reverse chronological)
% Please cite Meng et al.2017(at minimum) if this script contributes to 
% your publication.
%
% 1) UNDER REVIEW
%    H. L. White, W. Chen, N. M. Pugno, and S. Keten, "Rate-dependent size
%    effects govern the inverse thickness dependence of specific
%    penetration energy in nanoscale thin films"
% 2) H. L. White, A. Giuntoli, M. Fermen-Coker, and S. Keten, “Tailoring 
%    flake size and chemistry to improve impact resistance of graphene 
%    oxide thin films,” Carbon, vol. 215, p. 118382, Nov. 2023, 
%    doi: 10.1016/j.carbon.2023.118382.
% 3) T. Li, Z. Meng, and S. Keten, “Interfacial mechanics and viscoelastic 
%    properties of patchy graphene oxide reinforced nanocomposites,” 
%    Carbon, vol. 158, pp. 303–313, Mar. 2020, 
%    doi: 10.1016/j.carbon.2019.10.039.
% 4) Z. Meng et al., “A coarse-grained model for the mechanical behavior of 
%    graphene oxide,” Carbon, vol. 117, pp. 476–487, June 2017, 
%    doi: 10.1016/j.carbon.2017.02.061.
%
% Note: While dihedral terms are not used in the graphene oxide model,
% they are still identified in this script and can be output by 
% uncommenting lines 427, 436, and 543-548.
%
% MIT License
% 
% Copyright (c) 2025 Keten Research Group at Northwestern University
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% *************************** END DESCRIPTION *****************************

clear;
clc;

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BASIC INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assume units are Angstroms unless otherwise stated

% nx = number of repeat units in each flake's x (zigzag) direction
% nx = ceil[(zigzag_length+2.475)/2.477]
nx = 21; % Should be odd to avoid dangling beads

% ny = number of repeat units in each flake's y (armchair) direction
% ny = ceil[(armchair_length+2.864)/4.29]
ny = 12; % Should be even to avoid dangling beads

x_nsheetsperlayer = 3; % Number of flakes per layer in x (zigzag)
y_nsheetsperlayer = 3; % Number of flakes per layer in y (armchair)

nlayers = 6; % total number of layers

% Describe oxidation levels
pervac1 = 20; % percentage of beads B (hydroxyl oxidation)
pervac2 = 20; % percentage of beads C (epoxide oxidation)

% To randomly overlap flakes
% 0 = false, 1 = true
x_random_shift_status = 0;
y_random_shift_status = 0;

% To overlap flakes by a specific fraction
% NOTE: overridden if corresponding random_shift_status = 1
x_shift_factor = 0.5; 
y_shift_factor = 0.5;

% Distance between adjacent GO flakes
x_defectsize = 8;
y_defectsize = 8;

% Extra space at edges of simulation cell 
roomz = 4; roomx = x_defectsize/2; roomy = y_defectsize/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADVANCED INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
mapping_factor = 2; % mapping 4:1 --> 2; mapping 16:1 --> 4;
b0 = 1.43; % bond length graphene
b = b0*mapping_factor;
incrb = b*sin(pi/3); % ~2.4 A for 4:1 mapping
incry = 1.5*b; % ~ 4.2 A for 4:1 mapping

% Interlayer spacing based on nonbonded interactions (see publication 2
% for latest parameters)
% USE THIS FOR GO
sigma = 6.8; req = 2^(1/6)*sigma;
% USE THIS FOR G (no O)
%req = 3.46; sigma = req / (2^(1/6));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GO coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = []; y = []; z = [];

%Sample sheet for calculations of sheet width/length & box size
for j=1:ny
    for i=1:nx
        xprov = (i-1)*incrb;
        zprov = 0;
        if mod(j,2)==1 % j is odd
            if mod(i,2)==1 % i is odd
                yprov = b*cos(pi/3)+(j-1)*incry;
            else % i is even
                yprov = (j-1)*incry;
            end
        else % j is even
            if mod(i,2)==1 % odd
                yprov = b*(1+cos(pi/3))+(j-2)*incry;
            else
                yprov = b*(1+2*cos(pi/3))+(j-2)*incry;
            end
        end
        x = [x; xprov];
        y = [y; yprov];
        z = [z; zprov];
    end
end

sheet_length_A = max(x) - min(x);
sheet_width_A = max(y) - min(y);

minx_firstsheet = min(x);
maxx_firstsheet = max(x);
miny_firstsheet = min(y);
maxy_firstsheet = max(y);

% Erase sample sheet
x = []; y = []; z = [];

% Generate coordinates
for ysheet = 0:y_nsheetsperlayer-1
    for xsheet = 0:x_nsheetsperlayer-1
        for j=1:ny
            for i=1:nx

                xprov = (i-1)*incrb + xsheet*(sheet_length_A + x_defectsize);
                zprov = 0;

                if mod(j,2)==1 % j is odd
                    if mod(i,2)==1 % i is odd
                        yprov = b*cos(pi/3)+(j-1)*incry + ysheet*(sheet_width_A + y_defectsize);
                    else % i is even
                        yprov = (j-1)*incry + ysheet*(sheet_width_A + y_defectsize);
                    end

                else % j is even
                    if mod(i,2)==1 % odd
                        yprov = b*(1+cos(pi/3))+(j-2)*incry + ysheet*(sheet_width_A + y_defectsize);
                    else
                        yprov = b*(1+2*cos(pi/3))+(j-2)*incry + ysheet*(sheet_width_A + y_defectsize);
                    end
                end

                x = [x; xprov];
                y = [y; yprov];
                z = [z; zprov];
            end
        end
    end
end

go_coordinates = [x y z];

% Bonds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
go_bond = [];

for ysheet = 0:y_nsheetsperlayer-1
    for xsheet = 0:x_nsheetsperlayer-1
        % horizontal bonds
        for j=1:ny
            for i=1:nx-1
                go_bond = [go_bond; ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i+1];
            end
        end
        % vertical bonds
        for j=1:ny-1
            if mod(j,2)==1 % j is odd
                for i=1:2:nx
                    go_bond = [go_bond; ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i];
                end
            else
                for i=2:2:nx
                    go_bond = [go_bond; ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i];
                end
            end
        end
    end
end

% Angles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
go_angle = [];

for ysheet = 0:y_nsheetsperlayer-1
    for xsheet = 0:x_nsheetsperlayer-1
        % horizontal angles
        for j=1:ny
            for i=1:(nx-2)
                go_angle = [go_angle; ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i+1 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i+2];
            end
        end
        % vertical angles
        for j=1:ny-1
            if mod(j,2)==1 % j is odd
                for i=3:2:nx-1
                    go_angle = [go_angle;
                                ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i-1 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i; 
                                ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i+1 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i; 
                                ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i-1     ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i     ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i; 
                                ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i+1     ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i     ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i];
                end
            else % j is even
                for i=2:2:nx-2
                    go_angle = [go_angle; 
                                ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i-1 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i; 
                                ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i+1 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i; 
                                ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i-1     ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i     ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i; 
                                ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i+1     ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i     ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i];
                end
            end
        end
    end
end

% Dihedrals
go_dihedral = [];
for ysheet = 0:y_nsheetsperlayer-1
    for xsheet = 0:x_nsheetsperlayer-1
        % vertical dihedrals
        for j=1:2:ny-3 % type 1 odd
            for i=1:2:nx-2
                go_dihedral = [go_dihedral; ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i+1 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j+1)*nx+i+1];
            end
        end
        for j=2:2:ny-2 % type 1 even
            for i=2:2:nx-1
                go_dihedral = [go_dihedral; ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i+1 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j+1)*nx+i+1];
            end
        end

        for j=1:2:ny-1 % type 2 odd
            for i=1:2:nx-2
                go_dihedral = [go_dihedral; ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i+1 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i+1];
            end
        end
        for j=2:2:ny-2 % type 2 even
            for i=2:2:nx-1
                go_dihedral = [go_dihedral; ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i+1 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i+1];
            end
        end

        for j=1:2:ny-1 % type 3 odd
            for i=2:2:nx-1
                go_dihedral = [go_dihedral; ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i+1 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i+1 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i];
            end
        end
        for j=2:2:ny-2 % type 3 even
            for i=1:2:nx-2
                go_dihedral = [go_dihedral; ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i+1 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i+1 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i];
            end
        end

        for j=1:2:ny-1 % type 4 odd
            for i=2:2:nx-3
                go_dihedral = [go_dihedral; ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i+1 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i+1 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i+2];
            end
        end
        for j=2:2:ny-2 % type 4 even
            for i=1:2:nx-2
                go_dihedral = [go_dihedral; ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i+1 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i+1 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+j*nx+i+2];
            end
        end

        % horizontal dihedrals: type 5
        for j=1:ny
            for i=1:nx-3
                go_dihedral = [go_dihedral; ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i+1 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i+2 ysheet*x_nsheetsperlayer*nx*ny+xsheet*ny*nx+(j-1)*nx+i+3];
            end
        end

        go_bead_type = ones(length(go_coordinates),1);
        go_dihedral_type = ones(length(go_dihedral),1);
        go_angle_type = ones(length(go_angle),1);
        go_bond_type = ones(length(go_bond),1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% end of GO sheets %%%%%%%%%%%%%%%%%%%%%%%%%%


stack_coordinates = [];
for i=1:nlayers
    prov = [go_coordinates(:,1) go_coordinates(:,2) go_coordinates(:,3)+(i-1)*req];
    stack_coordinates = [stack_coordinates; prov];
end
comx = 0.5*(max(stack_coordinates(:,1))-min(stack_coordinates(:,1)));
comy = 0.5*(max(stack_coordinates(:,2))-min(stack_coordinates(:,2)));
stack_coordinates = [stack_coordinates(:,1)-comx stack_coordinates(:,2)-comy stack_coordinates(:,3)];

stack_bond = [];
stack_angle = [];
stack_dihedral = [];
for i=1:nlayers
    stack_bond = [stack_bond; go_bond+(i-1)*length(go_coordinates)];
    stack_angle = [stack_angle; go_angle+(i-1)*length(go_coordinates)];
    stack_dihedral = [stack_dihedral; go_dihedral+(i-1)*length(go_coordinates)];
end

stack_bead_type = repmat(go_bead_type,nlayers,1); % no dreiding type
stack_bond_type = repmat(go_bond_type,nlayers,1);
stack_angle_type = repmat(go_angle_type,nlayers,1);
stack_dihedral_type = repmat(go_dihedral_type,nlayers,1);

%%%%%%%%%%%%%%%%%%%%% introduce oxidized bead type %%%%%%%%%%%%%%%%%%%%%

rng('shuffle');
prov = randperm(nlayers*length(go_coordinates));
indecestypeb = prov(1:floor(pervac1/100*nlayers*length(go_coordinates))); %% pick up alchol beads
indecestypec = prov(floor(pervac1/100*nlayers*length(go_coordinates))+1:floor((pervac1+pervac2)/100*nlayers*length(go_coordinates))); %% pick up epoxide beads

indecesb = sort(indecestypeb);
indecesc = sort(indecestypec);

%%%%%%%%%%%%%%%%%%%%% change bead angle dihedral type %%%%%%%%%%%%%%%%%%%%%
stack_bead_type(indecesb) = [2];
stack_bead_type(indecesc) = [3];

bond_b_check = ismember(stack_bond,indecesb);
bond_c_check = ismember(stack_bond,indecesc);

stack_bond_type(any(bond_b_check,2)) = [2];% if either of the 2 beads in a bond is type B, then the bond type is 2
stack_bond_type(any(bond_c_check,2)) = [3]; % if either of the 2 beads in a bond is type C, then the bond type is 3

stack_angle_type(sum(ismember(stack_angle,indecesb),2)>1) = [2]; % if more than one bead in an angle is type B, then the angle type is 2
stack_angle_type(sum(ismember(stack_angle,indecesc),2)>1) = [3]; % if more than one bead in an angle is type B, then the angle type is 2

stack_dihedral_type(sum(ismember(stack_dihedral,indecesb),2)>2) = [2]; % if more than two beads in a dihedral is type B, then the dihedral type is 2
stack_dihedral_type(sum(ismember(stack_dihedral,indecesc),2)>2) = [3]; % if more than two beads in a dihedral is type B, then the dihedral type is 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MOLECULE TAGS
molecule_tag = [];
for i=1:nlayers*x_nsheetsperlayer*y_nsheetsperlayer
    flakesize = length(go_coordinates)/(x_nsheetsperlayer*y_nsheetsperlayer);
    molecule_tag = [molecule_tag; zeros(flakesize,1)+i];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% charge
charge = zeros(length(stack_coordinates),1);

% indexes
indexc = (1:length(stack_coordinates))';
indexb = (1:length(stack_bond))';
indexa = (1:length(stack_angle))';
indexd = (1:length(stack_dihedral))';

%Inputs for printing the data file
coordinates = [indexc molecule_tag stack_bead_type charge stack_coordinates];
bonds = [indexb stack_bond_type stack_bond];
angles = [indexa stack_angle_type stack_angle];
dihedrals = [indexd stack_dihedral_type stack_dihedral];

% inputs in this section are: coordinates, bonds, angles, dihedrals,
% impropers
% Write LAMMPS data file open the file with write permission
% inputs in this section are: coordinates, b, a,
% Write LAMMPS data file open the file with write permission

prov1 = num2str(nlayers);
al = num2str(pervac1);
ep = num2str(pervac2);
xspl = num2str(x_nsheetsperlayer);
yspl = num2str(y_nsheetsperlayer);
nx_str = num2str(nx);
ny_str = num2str(ny);

if x_random_shift_status == 1 || y_random_shift_status == 1

    if x_random_shift_status == 1 && y_random_shift_status == 1
        name = [ 'layers' prov1 '_' nx_str 'nx_' ny_str 'ny_' al 'a_' ep 'e_' xspl 'xspl_' yspl 'yspl_randxy.data'];
    elseif x_random_shift_status == 1 && y_random_shift_status ~= 1
        name = [ 'layers' prov1 '_' nx_str 'nx_' ny_str 'ny_' al 'a_' ep 'e_' xspl 'xspl_' yspl 'yspl_randx.data'];
    elseif x_random_shift_status ~= 1 && y_random_shift_status == 1
        name = [ 'layers' prov1 '_' nx_str 'nx_' ny_str 'ny_' al 'a_' ep 'e_' xspl 'xspl_' yspl 'yspl_randy.data'];
    end

else
    name = [ 'layers' prov1 '_' nx_str 'nx_' ny_str 'ny_' al 'a_' ep 'e_' xspl 'xspl_' yspl 'yspl.data'];
end

fid = fopen(name,'w');

fprintf(fid,'%s\n\n','Multilayered Graphene Oxide (Keten Research Group)');

natoms = length(coordinates);
nbonds = length(bonds);
nangles = length(angles);
ndihedrals = length(dihedrals);
nimpropers = 0;
fprintf(fid,'%d atoms\n',natoms);
fprintf(fid,'%d bonds\n',nbonds);
fprintf(fid,'%d angles\n\n',nangles);
% fprintf(fid,'%d dihedrals\n\n',ndihedrals);

natomtypes = 3;
nbondtypes = 3;
nangletypes = 3;
ndihedraltypes = 3;
fprintf(fid,'%d atom types\n',natomtypes);
fprintf(fid,'%d bond types\n',nbondtypes);
fprintf(fid,'%d angle types\n\n',nangletypes);
% fprintf(fid,'%d dihedral types\n\n',ndihedraltypes);

minx = min(stack_coordinates(:,1))-roomx;
maxx = max(stack_coordinates(:,1))+roomx;
miny = min(stack_coordinates(:,2))-roomy;
maxy = max(stack_coordinates(:,2))+roomy;
minz = min(stack_coordinates(:,3))-roomz;
maxz = max(stack_coordinates(:,3))+roomz;

num_indices = size(coordinates, 1);
num_indices_per_sheet = num_indices/nlayers;

% For predetermined shift amount
if x_shift_factor > 0
    xShiftBy = x_shift_factor*(maxx_firstsheet - minx_firstsheet)+x_defectsize/2;
else
    xShiftBy = 0;
end

if y_shift_factor > 0 
    yShiftBy = y_shift_factor*(maxy_firstsheet - miny_firstsheet)+y_defectsize/2;
else
    yShiftBy = 0;
end

if x_random_shift_status == 0 || y_random_shift_status == 0
    for whichsheet = 1:2:nlayers
        lowindex = num_indices_per_sheet*(whichsheet-1)+1;
        highindex = num_indices_per_sheet*whichsheet;
        for item = lowindex:highindex
            if x_random_shift_status == 0
                coordinates(item, 5) = coordinates(item, 5) + xShiftBy;
            end
            if y_random_shift_status == 0
                coordinates(item, 6) = coordinates(item, 6) + yShiftBy;
            end
        end
    end
end

% For random shift amount
if x_random_shift_status == 1 || y_random_shift_status == 1
    for whichsheet = 1:nlayers
        lowindex = num_indices_per_sheet*(whichsheet-1)+1;
        highindex = num_indices_per_sheet*whichsheet;

        if x_random_shift_status == 1
            x_shift_factor = round(rand(),2);
            xShiftBy = x_shift_factor*(maxx_firstsheet - minx_firstsheet)+x_defectsize/2;
        end

        if y_random_shift_status == 1
            y_shift_factor = round(rand(),2);
            yShiftBy = y_shift_factor*(maxy_firstsheet - miny_firstsheet)+y_defectsize/2;
        end

        for item = lowindex:highindex            
            if x_random_shift_status == 1
                coordinates(item, 5) = coordinates(item, 5) + xShiftBy;
            end
            if y_random_shift_status == 1
                coordinates(item, 6) = coordinates(item, 6) + yShiftBy;
            end
        end
    end
end

% the slack is correct for even nx and even ny
minx = min(stack_coordinates(:,1))-roomx;
maxx = max(stack_coordinates(:,1))+roomx;
miny = min(stack_coordinates(:,2))-roomy;
maxy = max(stack_coordinates(:,2))+roomy;
minz = min(stack_coordinates(:,3))-roomz;
maxz = max(stack_coordinates(:,3))+roomz;
fprintf(fid,'%.2f %.2f xlo xhi\n',minx,maxx);
fprintf(fid,'%.2f %.2f ylo yhi\n',miny,maxy);
fprintf(fid,'%.2f %.2f zlo zhi\n\n',minz,maxz);

fprintf(fid,'%s\n\n','Masses');
mass = [48; 65; 64]; % 1=sp2, 2=oxide, 3=filler
id = [1; 2; 3];
masses = [id mass];
for i=1:length(id)
    prov = masses(i,:);
    fprintf(fid,'%d %.2f\n',prov);
end

fprintf(fid,'\n%s\n\n','Atoms');
for i=1:length(coordinates)
    prov = coordinates(i,:);
    fprintf(fid,'%d %d %d %.2f %.2f %.2f %.2f\n',prov);
end

%print bonds
fprintf(fid,'\n%s\n\n','Bonds');
for i=1:length(bonds)
    prov = bonds(i,:);
    fprintf(fid,'%d %d %d %d\n',prov);
end

%print angles
fprintf(fid,'\n%s\n\n','Angles');
for i=1:length(angles)
    prov = angles(i,:);
    fprintf(fid,'%d %d %d %d %d\n',prov);
end

%print dihedrals
% fprintf(fid,'\n%s\n\n','Dihedrals');
% for i=1:length(dihedrals)
%     prov = dihedrals(i,:);
%     fprintf(fid,'%d %d %d %d %d %d\n',prov);
% end

fclose(fid);

toc;
