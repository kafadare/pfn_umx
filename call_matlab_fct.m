function out = call_matlab_fct(index, file)

cd '/cbica/projects/bgd-pfn' ;

%dataline = [index+1 index+1];

%opts = detectImportOptions(file)
%opts.DataLines = dataline;

%A = csvread(file, index,0, [index 0 index 1]);
T = readtable(file);
%n = A(1);
%k = A(2);
n = T.Vertex_No;
k = T.PFN;
out = get_vertex_col(n,k);
end
