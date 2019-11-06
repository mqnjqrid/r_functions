# rowcal is the box at the intersection of rownames and column names
mean_var_table = function(table1, table2, rowcol = "", colname = colnames(table1), rowname = rownames(table1)) {
  if(nrow(table1) != nrow(table2) | ncol(table1) != ncol(table2)){
    print("Error: table dimension do not match")
    stop()
    }
  for(rownum in 0:nrow(table1)) {
    if(rownum == 0){
      cat("\\begin{tabular}{|c|"}, rep('r', ncol(table1)), "\\hline\n", sep = '')
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
