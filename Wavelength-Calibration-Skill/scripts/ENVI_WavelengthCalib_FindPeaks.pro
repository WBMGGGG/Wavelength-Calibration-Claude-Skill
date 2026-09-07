; ENVI 波长标定批量找峰脚本
; 运行方式：ENVI→文件→打开脚本→运行本脚本
PRO WavelengthCalib_FindPeaks
  COMPILE_OPT idl2
  ; ========= 用户配置 =========
  dataRoot = ''  ; 日期根文件夹路径
  approx_a = 0.4
  approx_b = 400.0
  lambda_std = [404.66,435.83,546.07,576.96,579.07,585.25,594.48,614.31,$
    696.54,703.24,724.52,743.89,750.39,763.51,772.42,794.82,$
    800.62,811.53,826.45,842.46,852.14,912.30]
  lampRow = {hg:[0,1,2], ne:[3,4,5,6,7,9,10,11], ar:[8,12,13,14,15,16,17,18,19,20,21]}
  ; ========================

  ENVI, /RESTORE_BASE
  cd, dataRoot
  subDirs = FILE_SEARCH('*', /TEST_DIRECTORY, /NOSORT)
  outTable = LIST()

  FOREACH dir, subDirs DO BEGIN
    IF ~FILE_TEST(dir+'.hdr') THEN CONTINUE
    lamp = ''
    IF STRPOS(dir,'HG') GE 0 THEN lamp='hg'
    IF STRPOS(dir,'NE') GE 0 THEN lamp='ne'
    IF STRPOS(dir,'AR') GE 0 THEN lamp='ar'
    IF lamp EQ '' THEN CONTINUE

    ENVI_OPEN_FILE, dir+'.hdr', R_FID=fid
    ENVI_QUERY, fid, NSAMPLES=samples, NLINES=lines, DATA_TYPE=dt
    data = ENVI_GET_DATA(fid)
    centerLine = FIX(lines/2)
    spec = data[*, centerLine]

    rows = lampRow.(lamp)
    pixelRes = MAKE_ARRAY(N_ELEMENTS(rows), /FLOAT, VALUE=!VALUES.F_NAN)

    FOR k=0, N_ELEMENTS(rows)-1 DO BEGIN
      tl = lambda_std[rows[k]]
      cPix = (tl - approx_b) / approx_a
      sr = FIX(2 / approx_a)
      st = 0 > (FIX(cPix)-sr) ? 0 : (FIX(cPix)-sr)
      ed = (samples-1) < (FIX(cPix)+sr) ? (samples-1) : (FIX(cPix)+sr)
      seg = spec[st:ed]

      peaks = PEAKFIND(seg, MINPEAKHEIGHT=MAX(seg)*0.3)
      IF N_ELEMENTS(peaks) EQ 0 THEN CONTINUE

      minDiff = 1e9
      selPix = -1
      FOREACH p, peaks DO BEGIN
        pp = p + st
        lamEst = approx_a * pp + approx_b
        diff = ABS(lamEst - tl)
        IF diff LE 2.0 AND diff LT minDiff THEN BEGIN
          minDiff = diff
          selPix = pp
        ENDIF
      ENDFOREACH

      IF selPix GE 0 THEN BEGIN
        IF selPix GE 2 AND selPix LE samples-3 THEN BEGIN
          winTop = spec[selPix-1:selPix+1]
          dMax = MAX(winTop) - MIN(winTop)
          IF dMax LT 0.02*MAX(spec) THEN BEGIN
            x = [-1,0,1]
            y = winTop
            fit = GAUSSFIT(x, y)
            selPix = selPix + fit[1]
          ENDIF
        ENDIF
        pixelRes[k] = selPix
      ENDIF
    ENDFOR

    outTable.Add, {dir:dir, lamp:lamp, pixels:pixelRes}
    ENVI_CLOSE, fid
  ENDFOREACH

  ; 输出结果到控制台，可复制回填Excel
  PRINT, '===== 峰检测结果 ====='
  FOREACH item, outTable DO BEGIN
    PRINT, '文件夹: ', item.dir
    PRINT, '像素: ', item.pixels
  ENDFOREACH
END
