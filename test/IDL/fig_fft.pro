!p.charsize=1.8
!p.thick=5
!x.thick=5
!y.thick=5

log_spacex=9
log_spacey=2

set_device,'fft_G28.eps',/port

;filename='../G8/run6[68]*/GM/cut.outs'
filename='../G28/run_G28_PIC_0810/GM/box*.outs'
;filename='../G28/run_PIC_G28_FaceInnerBC_599s_0603/GM/box*.outs'
colors=[255,250,100]
smooths=[5,5,5,5]
extract_galileo_times,firstpict=60, timeshift0=201, stretch=1.0, $
                      freqrange=[0,0.3],flyby=8

close_device,/pdf

set_default_values

;;; set_device,'fig_fft_inside.eps',/port
;;;
;;; extract_galileo_times,filename='output/run6[68]*/GM/cut.outs', $
;;;                       firstpict=60, timeshift0=261, stretch=1.06, $
;;;                       xrange=[952,962], $
;;;                       brange=[1e-6,1e4], $
;;;                       freqrange=[0,0.3], $
;;;                       colors=[255,250,100]
;;;
;;; close_device,/pdf
