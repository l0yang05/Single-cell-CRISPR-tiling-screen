## functions
add_category <-function(feature_call, anno.tbl){
  ## if feature_call does not exisit in the anno,tbl, add feature_call as a new category
  idx = match(feature_call, anno.tbl$id)
  if(is.na(idx)){
    cat.info<-feature_call
  }
  else{cat.info <- anno.tbl$category[idx]}
  cat.info
}
