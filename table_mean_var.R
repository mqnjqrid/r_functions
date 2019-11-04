# rowcal is the box at the intersection of rownames and column names
mean_var_table = function(table1, table2, rowcol = "", colname = colnames(table1), rowname = rownames(table1)) {
for(rownum in 0:nrow(table1)) {
  if(rownum == 0){
    cat("\\begin{tabular}{|c|rrrr|} \\hline\n")
    cat(rowcol, " & ")
    cat(colname, sep = " & ")
    cat("\\\\ \\hline\n")
  }else{
    if(rownum == nrow(table1)){
      cat("\\hline\n")
    }
    cat(rowname[rownum])
    for(colnum in 1:ncol(table1)){
       cat(" & ", table1[rownum, colnum], ' (', table2[rownum, colnum], ')', sep = '')
    }
    cat("\\\\ \n")
  }
  if(rownum == nrow(table1)){
    cat("\\hline\n \\end{tabular}")
  }
}
}
