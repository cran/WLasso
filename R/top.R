top <-
function(x, thresh, sorted_vect){
  v=sorted_vect[thresh]
  ifelse(abs(x)>v,x,0)
}
