;Generate m_tilda

ns1= random_uniform(-1, 1, 500000) 
ns2= random_uniform(-1, 1, 500000) 

;print(ns)

m_til1=new((/dimsizes(ns1)/),float)
m_til2=m_til1

;intial values

m_til1(0)=-1.
m_til1(1)=0.5

m_til2(0)=1.
m_til2(1)=0.5

  do t=2, dimsizes(ns1)-1
    m_til1(t)=0.6*m_til1(t-1) - 0.3*m_til1(t-2) + ns1(t)
    m_til2(t)=0.6*m_til2(t-1) - 0.3*m_til2(t-2) + ns2(t)
  end do

;m_til1=random_uniform(-1, 1, 30000)
;m_til2=random_uniform(-1, 1, 30000)


;generate synthetic data u and m

ut1=new((/dimsizes(ns1)/),float)
u1 =new((/dimsizes(ns1)/),float)
m1 =new((/dimsizes(ns1)/),float)

ut2=new((/dimsizes(ns1)/),float)
u2 =new((/dimsizes(ns1)/),float)
m2 =new((/dimsizes(ns1)/),float)

F1=1./11.
F2=1./11.

;imposed feedbacks

b11 =  0.04 
b12 =  0.06  *0
b21 = -0.025  *0
b22 =  0.00

;real

;b11 = 0.03775075 
;b12 = 0.05863772   *0
;b21 = -0.01757893  *0
;b22 =  0.0264


u1(0)= 10.   ;initial value
u2(0)= 10.   ;initial value

  do t=1, dimsizes(ns1)-1
    m1 (t) = m_til1(t) + b11*u1(t-1) + b12*u2(t-1)
    m2 (t) = m_til2(t) + b21*u1(t-1) + b22*u2(t-1)
    ut1(t) = m1(t) - F1*u1(t-1)
    ut2(t) = m2(t) - F2*u2(t-1)
    u1 (t) = ut1(t)+u1(t-1)
    u2 (t) = ut2(t)+u2(t-1)        
  end do



;cross correlation

  mxlag    = 60
  lag      = ispan(-mxlag, mxlag, 1)

;u1u1

  xy2   = esccr(u1,u1,mxlag)
  yx2   = esccr(u1,u1,mxlag)
  ccr2  = new ((/2*mxlag+1/), float)
  ccr2(0:mxlag-1) = yx2(1:mxlag:-1)
  ccr2(mxlag:)    = xy2(0:mxlag)
  ccr2!0="lag"
  ccr2&lag=lag
  printVarSummary(ccr2)

;u2u2

  xy3   = esccr(u2,u2,mxlag)
  yx3   = esccr(u2,u2,mxlag)
  ccr3  = new ((/2*mxlag+1/), float)
  ccr3(0:mxlag-1) = yx3(1:mxlag:-1)
  ccr3(mxlag:)    = xy3(0:mxlag)
  ccr3!0="lag"
  ccr3&lag=lag
  printVarSummary(ccr3)

 ccr23  = new ((/2,2*mxlag+1/), float)
 ccr23(0,:)=ccr2
 ccr23(1,:)=ccr3

;u1u2
  xy4   = esccr(u1,u2,mxlag)
  yx4   = esccr(u2,u1,mxlag)
  ccr4  = new ((/2*mxlag+1/), float)
  ccr4(0:mxlag-1) = yx4(1:mxlag:-1)
  ccr4(mxlag:)    = xy4(0:mxlag)
  ccr4!0="lag"
  ccr4&lag=lag
  printVarSummary(ccr4)
  ccr4=ccr4*-1

;m1m2
  xy5   = esccr(m2,m1,mxlag)
  yx5   = esccr(m1,m2,mxlag)
  ccr5  = new ((/2*mxlag+1/), float)
  ccr5(0:mxlag-1) = yx5(1:mxlag:-1)
  ccr5(mxlag:)    = xy5(0:mxlag)
  ccr5!0="lag"
  ccr5&lag=lag
  printVarSummary(ccr5)

;m1z1
  xy6   = esccr(u1,m1,mxlag)
  yx6   = esccr(m1,u1,mxlag)
  ccr6  = new ((/2*mxlag+1/), float)
  ccr6(0:mxlag-1) = yx6(1:mxlag:-1)
  ccr6(mxlag:)    = xy6(0:mxlag)
  ccr6!0="lag"
  ccr6&lag=lag
  printVarSummary(ccr6)

;m2z2
  xy7   = esccr(u2,m2,mxlag)
  yx7   = esccr(m2,u2,mxlag)
  ccr7  = new ((/2*mxlag+1/), float)
  ccr7(0:mxlag-1) = yx7(1:mxlag:-1)
  ccr7(mxlag:)    = xy7(0:mxlag)
  ccr7!0="lag"
  ccr7&lag=lag
  printVarSummary(ccr7)

;m1z2
  xy8   = esccr(u2,m1,mxlag)
  yx8   = esccr(m1,u2,mxlag)
  ccr8  = new ((/2*mxlag+1/), float)
  ccr8(0:mxlag-1) = yx8(1:mxlag:-1)
  ccr8(mxlag:)    = xy8(0:mxlag)
  ccr8!0="lag"
  ccr8&lag=lag
  printVarSummary(ccr8)


;m2z1
  xy9   = esccr(u1,m2,mxlag)
  yx9   = esccr(m2,u1,mxlag)
  ccr9  = new ((/2*mxlag+1/), float)
  ccr9(0:mxlag-1) = yx9(1:mxlag:-1)
  ccr9(mxlag:)    = xy9(0:mxlag)
  ccr9!0="lag"
  ccr9&lag=lag
  printVarSummary(ccr9)


print(u1(0))
print(u1(1))

print(u2(0))
print(u2(1))

time=ispan(1,dimsizes(u1),1)
dtu1=center_finite_diff_n (u1, time, False, 0, 0) 
print(dtu1(0))

dtu2=center_finite_diff_n (u2, time, False, 0, 0) 
print(dtu2(0))



;;;;;;;; 



;;; Autocor


acr_z1_pos = esacr(u1, mxlag)
acr_z1_neg = esacr(u1, mxlag)
acr_z1 = new(2*mxlag+1, double)
acr_z1(0:mxlag-1) = acr_z1_neg(1:mxlag:-1)
acr_z1(mxlag:) = acr_z1_pos

acr_z2_pos = esacr(u2, mxlag)
acr_z2_neg = esacr(u2, mxlag)
acr_z2 = new(2*mxlag+1, double)
acr_z2(0:mxlag-1) = acr_z2_neg(1:mxlag:-1)
acr_z2(mxlag:) = acr_z2_pos

acr_m1_pos = esacr(m1, mxlag)
acr_m1_neg = esacr(m1, mxlag)
acr_m1 = new(2*mxlag+1, double)
acr_m1(0:mxlag-1) = acr_m1_neg(1:mxlag:-1)
acr_m1(mxlag:) = acr_m1_pos

acr_m2_pos = esacr(m2, mxlag)
acr_m2_neg = esacr(m2, mxlag)
acr_m2 = new(2*mxlag+1, double)
acr_m2(0:mxlag-1) = acr_m2_neg(1:mxlag:-1)
acr_m2(mxlag:) = acr_m2_pos

sum_acr_z1_z2 = sum(acr_z1 * acr_z2)
sum_acr_m1_z1 = sum(acr_m1 * acr_z1)
sum_acr_m2_z2 = sum(acr_m2 * acr_z2)
sum_acr_m1_z2 = sum(acr_m1 * acr_z2)
sum_acr_m2_z1 = sum(acr_m2 * acr_z1)

var_z1_z2 = new((2*mxlag+1), double)
do i = -mxlag, mxlag
	var_z1_z2(i + mxlag) = (1.0 / (dimsizes(u1) - abs(i))) * sum_acr_z1_z2
end do
s_z1_z2 = sqrt(var_z1_z2)

var_m1_z1 = new((2*mxlag+1), double)
do i = -mxlag, mxlag
	var_m1_z1(i + mxlag) = (1.0 / (dimsizes(u1) - abs(i))) * sum_acr_m1_z1
end do
s_m1_z1 = sqrt(var_m1_z1)

var_m2_z2 = new((2*mxlag+1), double)
do i = -mxlag, mxlag
	var_m2_z2(i + mxlag) = (1.0 / (dimsizes(u1) - abs(i))) * sum_acr_m2_z2
end do
s_m2_z2 = sqrt(var_m2_z2)

var_m1_z2 = new((2*mxlag+1), double)
do i = -mxlag, mxlag
	var_m1_z2(i + mxlag) = (1.0 / (dimsizes(u1) - abs(i))) * sum_acr_m1_z2
end do
s_m1_z2 = sqrt(var_m1_z2)

var_m2_z1 = new((2*mxlag+1), double)
do i = -mxlag, mxlag
	var_m2_z1(i + mxlag) = (1.0 / (dimsizes(u1) - abs(i))) * sum_acr_m2_z1
end do
s_m2_z1 = sqrt(var_m2_z1)



X = ispan(-mxlag, mxlag, 1)

xp_z1_z2 = new(2 * (2*mxlag+1), double)
yp_z1_z2 = new(2 * (2*mxlag+1), double)

do i = 0, 2*mxlag
	yp_z1_z2(i) = s_z1_z2(i)
	yp_z1_z2(2 * (2*mxlag+1) - 1 - i) = -s_z1_z2(i)
	xp_z1_z2(i) = X(i)
	xp_z1_z2(2 * (2*mxlag+1) - 1 - i) = X(i)
end do

xp_m1_z1 = new(2 * (2*mxlag+1), double)
yp_m1_z1 = new(2 * (2*mxlag+1), double)

do i = 0, 2*mxlag
	yp_m1_z1(i) = s_m1_z1(i)
	yp_m1_z1(2 * (2*mxlag+1) - 1 - i) = -s_m1_z1(i)
	xp_m1_z1(i) = X(i)
	xp_m1_z1(2 * (2*mxlag+1) - 1 - i) = X(i)
end do

xp_m2_z2 = new(2 * (2*mxlag+1), double)
yp_m2_z2 = new(2 * (2*mxlag+1), double)

do i = 0, 2*mxlag
	yp_m2_z2(i) = s_m2_z2(i)
	yp_m2_z2(2 * (2*mxlag+1) - 1 - i) = -s_m2_z2(i)
	xp_m2_z2(i) = X(i)
	xp_m2_z2(2 * (2*mxlag+1) - 1 - i) = X(i)
end do

xp_m1_z2 = new(2 * (2*mxlag+1), double)
yp_m1_z2 = new(2 * (2*mxlag+1), double)

do i = 0, 2*mxlag
	yp_m1_z2(i) = s_m1_z2(i)
	yp_m1_z2(2 * (2*mxlag+1) - 1 - i) = -s_m1_z2(i)
	xp_m1_z2(i) = X(i)
	xp_m1_z2(2 * (2*mxlag+1) - 1 - i) = X(i)
end do

xp_m2_z1 = new(2 * (2*mxlag+1), double)
yp_m2_z1 = new(2 * (2*mxlag+1), double)

do i = 0, 2*mxlag
	yp_m2_z1(i) = s_m2_z1(i)
	yp_m2_z1(2 * (2*mxlag+1) - 1 - i) = -s_m2_z1(i)
	xp_m2_z1(i) = X(i)
	xp_m2_z1(2 * (2*mxlag+1) - 1 - i) = X(i)
end do


printVarSummary(xp_m2_z1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 wks   = gsn_open_wks ("x11","cross_feedback_nonpropagating")              
 plot = new(6,graphic)                          ; create a plot array
  
  res          = True
  res@gsnDraw  = False                          ; don't draw
  res@gsnFrame = False   
  res@vpHeightF = 0.2        ; Changes the aspect ratio
  res@vpWidthF  = 0.32    
  res@gsnYRefLine           = 0.0 
  res@gsnXRefLine           = 0.0 
  res@gsnStringFontHeightF  = 0.02  

  res@tmXBLabelFontHeightF = 0.015                 ; resize tick labels
  res@tmYLLabelFontHeightF = 0.015  

 res@xyDashPattern  = 0  
 res@xyLineThicknesses = (/  2.0,   2.0/)          ; make second line thicker
 res@xyLineColors      = (/"blue","red"/)          ; change line color

 res@tiXAxisString      = "Lag (days)"      ; add an axis title
 res@tiXAxisFontHeightF=0.015
 res@tiYAxisString      = " "      ; add an axis title

res@trYMinF=-0.2
res@trYMaxF=1.0

res@tmXBMajorOutwardLengthF = 0.0
  res@tmYLMajorOutwardLengthF = 0.0
  res@tmXBMinorOutwardLengthF = 0.0
  res@tmYLMinorOutwardLengthF = 0.0

  ;res@gsnYRefLine           = (/ 0.001 /)
  ;res@gsnAboveYRefLineColor = (/ "gray90"/)
  ;res@gsnBelowYRefLineColor = (/ "gray90"/)

 res@gsnCenterString     = "(a) z~B~1~N~z~B~1~N~, z~B~2~N~z~B~2~N~ "
 plot(0)  = gsn_csm_xy (wks,lag,ccr23,res) 


res@trYMinF = -0.2
res@trYMaxF =  0.2

 res@xyLineColors      = (/"black","red"/)  

 res@gsnCenterString     = "(b) z~B~2~N~z~B~1~N~"
 plot(1)  = gsn_csm_xy (wks,lag,ccr4,res)

res@trYMaxF=0.55


res@gsnCenterString     = "(c) m~B~1~N~z~B~1~N~ "
;res@tiYAxisString      = "b (day~S4~-1~N~)"      ; add an axis title
 plot(2)  = gsn_csm_xy (wks,lag,ccr6,res)


res@gsnCenterString     = "(d) m~B~2~N~z~B~2~N~ "
;res@tiYAxisString      = "b (day~S4~-1~N~)"      ; add an axis title
 plot(3)  = gsn_csm_xy (wks,lag,ccr7,res)

res@trYMaxF =  0.2

res@gsnCenterString     = "(e) m~B~1~N~z~B~2~N~"
;res@tiYAxisString      = "b (day~S4~-1~N~)"      ; add an axis title
 plot(4)  = gsn_csm_xy (wks,lag,ccr8,res)

res@gsnCenterString     = "(f) m~B~2~N~z~B~1~N~ "
;res@tiYAxisString      = "b (day~S4~-1~N~)"      ; add an axis title
plot(5)  = gsn_csm_xy (wks,lag,ccr9,res)

 res@tiYAxisString      = "Corr. Coefficient"      ; add an axis title
 res@gsnCenterString     = "(g) m~B~1~N~m~B~2~N~ "
 ;plot(6)  = gsn_csm_xy (wks,lag,ccr5,res)


gsres = True
gsres@gsFillColor = "gray60"
gsres@gsFillOpacityF = 0.25

;dummy = gsn_add_polygon(wks,plot,xp,yp,gsres)

dummy_z1_z2 = gsn_add_polygon(wks,plot(1),xp_z1_z2*3,yp_z1_z2*3,gsres)
dummy_m1_z1 = gsn_add_polygon(wks,plot(2),xp_m1_z1*3,yp_m1_z1*3,gsres)
dummy_m2_z2 = gsn_add_polygon(wks,plot(3),xp_m2_z2*5,yp_m2_z2*5,gsres)
dummy_m1_z2 = gsn_add_polygon(wks,plot(4),xp_m1_z2*3,yp_m1_z2*3,gsres)
dummy_m2_z1 = gsn_add_polygon(wks,plot(5),xp_m2_z1*3,yp_m2_z1*3,gsres)

  txres                         = True
  txres@txPerimOn             = False
  txres@txFontHeightF         = 0.015
  txres@txBackgroundFillColor = "White"             
  gsn_text_ndc(wks,"z~B~1~N~ leads",0.80,.76,txres)  
  gsn_text_ndc(wks,"z~B~1~N~ leads",0.40,.42,txres)  
  gsn_text_ndc(wks,"z~B~2~N~ leads",0.80,.42,txres) 
  gsn_text_ndc(wks,"z~B~2~N~ leads",0.40,.12,txres)  
  gsn_text_ndc(wks,"z~B~1~N~ leads",0.80,.12,txres)

  resP                    = True                 ; modify the panel plot
  resP@gsnPanelMainString = " "     ; set main title
  resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelXWhiteSpacePercent = 5
  resP@gsnPanelRowSpec = True 
  resP@gsnPanelCenter  = False
  gsn_panel(wks,plot,(/2,2,2/),resP)  
