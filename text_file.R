
txt_lines <- function (fName, col_names = NULL,...){
 vecs <- list(...)
 tmp <- do.call(expand.grid, vecs)
 if(!is.null(col_names)){
   write.table(tmp, file = paste0(fName, ".csv"), row.names = F, col.names = col_names, sep=",")
 } else {
  write.table(tmp, file = paste0(fName, ".csv"), row.names = F, col.names = F, sep=",")
 }
 return(tmp)
}
 
# lines <- matrix(NA,59412,17)
# lines <- matrix(NA,59412,17)
# 
# 
# for (n in 1:59412) {
#   for (k in 1:17){
#     lines[n,k] = paste0("qsub matlab_test.sh \"" ,n, "\" \"" ,k, "\"" )
#   
#   }
# }
# 
# all = paste(lines, sep = "\n")
# fcon <- file("matlab_columns.txt")
# write(lines, fcon)
# close(fcon)
# 

#paste0("matlab -r 'try get_vertex_col(" ,n, "," ,k, "); catch; end; quit' " , collapse = "")


