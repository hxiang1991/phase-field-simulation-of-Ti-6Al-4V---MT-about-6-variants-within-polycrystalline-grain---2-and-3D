function [  ]= write_vtk_grid_values(nx, ny, sym_cor_mat, istep, data1, data2, data3, data4, data5)         

format short;

%% == open output file

fname = sprintf('time_%d_with_elasticity.vtk', istep);
out = fopen(fname, 'w');

npoint = nx* ny;

%% == reshape the Nx-Ny-Nz matrix into Nx* Ny* Nz-1 array for convenient output   

data1 = reshape(data1, nx* ny, 1); data2 = reshape(data2, nx* ny, 1);
data3 = reshape(data3, nx* ny, 1); data4 = reshape(data4, nx* ny, 1);
data5 = reshape(data5, nx* ny, 1);

%% == start writing ASCII VTK file:

% == header of VTK file

fprintf(out, '# vtk DataFile Version 2.0\n');
fprintf(out, 'time_10.vtk\n');
fprintf(out, 'ASCII\n');
fprintf(out, 'DATASET STRUCTURED_GRID\n');

% == coords of grid points:

fprintf(out, 'DIMENSIONS %5d  %5d  %5d\n', nx, ny, 1);
fprintf(out, 'POINTS %7d   float\n', npoint);

Array_sym_cor_mat = reshape(sym_cor_mat.', [], 1);

formatSpec = '%14.6e  %14.6e  %14.6e\n';
fprintf(out, formatSpec, Array_sym_cor_mat);

% == write grid point values:

fprintf(out,'POINT_DATA %5d\n',npoint);

fprintf(out,'SCALARS Variants  float  1\n');

fprintf(out,'LOOKUP_TABLE default\n');

fprintf(out,'%14.6e\n', data1);                % output grid value of Variants     

fprintf(out,'SCALARS Elastic_energy  float  1\n');

fprintf(out,'LOOKUP_TABLE default\n');

fprintf(out,'%14.6e\n', data2);                % output grid value of elastic energy density   

fprintf(out,'SCALARS s11  float  1\n');

fprintf(out,'LOOKUP_TABLE default\n');

fprintf(out,'%14.6e\n', data3);                % output grid value of s11, unit: MPa  

fprintf(out,'SCALARS s22  float  1\n');

fprintf(out,'LOOKUP_TABLE default\n');

fprintf(out,'%14.6e\n', data4);                % output grid value of s22, unit: MPa 

fprintf(out,'SCALARS s12  float  1\n');

fprintf(out,'LOOKUP_TABLE default\n');

fprintf(out,'%14.6e\n', data5);                % output grid value of s12, unit: MPa  

fclose(out);

end %endfunction

      
