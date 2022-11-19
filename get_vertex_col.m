

function out_T = get_vertex_col(n,k)

cd '/cbica/projects/bgd-pfn' ;

twins_T = readtable('/cbica/projects/bgd-pfn/ABCD_T_pairsample_long.csv');

path = strcat(pwd(), '/subj_PFNs');

dir_subj = dir(path);

file_ext = '/IndividualParcel_Final_sbj1_comp17_alphaS21_1_alphaL300_vxInfo1_ard0_eta0/final_UV.mat';

fName = sprintf('PFN%dV%d.csv', n, k);

colname = sprintf('PFN%dV%d', n, k);

out = cell(height(twins_T), 2);

for i = 1:height(twins_T)
subj = char(twins_T.subjectkey(i));
subj_no = extractAfter(subj,"NDAR_INV");
filepath = strcat(path, '/sub-NDARINV', subj_no, file_ext);
struct = load(filepath);
out{i,1} = subj;
out{i,2} = struct.V{1}(k,n); % first column vertex, second column PFN no
end

cd  '/cbica/projects/bgd-pfn/vertex_columns';

out_T = cell2table(out,"VariableNames",{'subjectkey' colname});
writetable(out_T,fName);
end
