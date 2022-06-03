#' Approximate within points
#'
#' @description Calculation of all points in between a given start point and
#' end point for a given total number of points
#'
#' @param x.B starting point (begin)
#' @param x.E end point
#' @param time.points number of points on total
#'
#' @return all exact points
#' @export
#'
calc.X.values<-function(x.B, x.E, time.points){

  i<-0:(time.points-1)
  x<-x.B+(x.E-x.B)/(time.points-1)*i
  return(x)
}




#' Find approximate fidelity with a logarithmic function
#'
#' @description Calculates for given start Fidelity and end Fidelity the within
#' Fidelity for a given number of time points using a logarithmic function to approximate Fidelity in between,
#' where with the par.slope stretching and compression of the function can be
#' influenced. par.slope has to be greater than 0, near zero then the values are
#' almost constant, the greater then the values are going against the linear line.
#'
#' @param time.points number of time points
#' @param Fid.End Fidelity at the end (final time point)
#' @param Fid.T1 Fidelity at the begin (starting time point)
#' @param par.slope Slope for the logarithmic function has to be greater 0, the greater towards 1 the more linear is the function
#'
#' @return a matrix: first column with time points, second column corresponding Fidelity for the time points
#' @export
#'
find.Fidelity.log<-function(time.points, Fid.End, Fid.T1, par.slope=1){

  #using Logarithmusfunktion, Basis ist e

  i<-0:(time.points-1)
  ##Reference Logarithmus-Funktion
  #Start-X-Wert
  ref.x1<-0.005#?ndert sich in allen log. Kurven
  #End-X-Wert
  ref.x2<-1#bleibt bei allen Kurven (ist deren Schnittpunkt)
  #merke Referenz-Y-Startwert f?r andere Funktionen
  ref.y1<-log(ref.x1)

  #Streckungs/Stauchungsparameter der log Funktions
  m<-par.slope
  #Berechnung neuen Start-X-Wert, hat den selben Funktionswert wie in der Referenzfunktion
  x.B.neu<-(exp(1)^(ref.y1))^(1/m)
  #Berechnung f?r alle Zeitpunkte x-Werte
  x<-calc.X.values(x.B=x.B.neu, x.E=ref.x2, time.points)
  #Berechne Funktionswerte der Funktion zu x-Werten
  f<-m*log(x)
  #NormierungsFaktor f?r Funktionswerte zum Mappen auf [Fid_B,Fid_E]
  norm<-(Fid.End-Fid.T1)/(m*log(x[time.points])-m*log(x[1]))
  #Normierung der Funktionswerte zum Mappen in neuen Funktionsbereich
  f.norm<-Fid.End+norm*f
  res<-cbind(time=1:time.points, Fidelity.Prozent=round(f.norm*100))
  return(res)
}



#' Find approximate fidelity with a linear function
#'
#' @description calculates for given start Fidelity and end Fidelity the within
#' Fidelity for a given number of time points using a linear function
#'
#' @param time.points number of time points
#' @param Fid.End Fidelity at the end (final time point)
#' @param Fid.T1 Fidelity at the begin (starting time point)
#'
#' @return a matrix: first column with time points, second column corresponding Fidelity for the time points
#' @export
#'
find.Fidelity.linear<-function(time.points, Fid.End, Fid.T1){

  m<-(Fid.T1-Fid.End)/(1-time.points)
  n<-(Fid.End-time.points*Fid.T1)/(1-time.points)
  x<-1:time.points

  f<-m*x+n
  res<-cbind(time=1:time.points, Fidelity.Prozent=round((m*x+n)*100)
  )
  return(res)
}

#' Find approximate fidelity with an exponential function
#'
#' @description Calculates for given start Fidelity and end Fidelity the within Fidelity
#' for a given number of time points using an exponential function,
#' where with the par.slope stretching and compression of the function can be influenced.
#' Here the exponential function is derived as the reflection of the logarithm function on the straight line
#' par.slope has to be greater than 0, near zero then the values are almost constant, the greater then the values are going against the linear line
#'
#' @param time.points number of time points
#' @param Fid.End Fidelity at the end (final time point)
#' @param Fid.T1 Fidelity at the begin (starting time point)
#' @param par.slope Slope for the logarithmic function has to be greater 0, the greater towards 1 the more linear is the function
#'
#' @return a matrix: first column with time points, second column corresponding Fidelity for the time points
#' @export
#'
find.Fidelity.exp<-function(time.points, Fid.End, Fid.T1, par.slope=1){

  #Berechnung Logarithmusfunktion
  res.log<-find.Fidelity.log(time.points, Fid.End, Fid.T1, par.slope=par.slope)
  #Berechnung Spiegelgeraden
  res.linear<-find.Fidelity.linear(time.points, Fid.End, Fid.T1)
  #Differenzberechnung
  diff.exp<-rev(res.linear[1,"Fidelity.Prozent"]+
                  (res.linear[time.points,"Fidelity.Prozent"]-res.log[,"Fidelity.Prozent"]))
  res<-cbind(time=1:time.points,
             Fidelity.Prozent= diff.exp)
  return(res)
}
