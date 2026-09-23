#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
orchestrator_clonodynamics.py
=============================

Single-file interactive master orchestrator for the ClonoDynamics project.

WHAT IS CLONODYNAMICS?
----------------------
ClonoDynamics is a modular framework for studying longitudinal T-cell receptor
(TCR) repertoire dynamics from paired technical replicates. The workflow keeps
separate the biological longitudinal signal, measurement/observation noise,
pseudo-longitudinal technical-reference behavior, synthetic positive controls,
and downstream robustness/sensitivity analyses.

The production architecture is organized into five analysis blocks:

    01_core
        Genuine longitudinal data, Steps 1-6.

    02_pseudo_reference
        Pseudo-longitudinal technical reference, including Steps 7-8.

    03_validation
        Genuine longitudinal forward/fluctuation/temporal analyses, Steps 9-13,
        plus the support-conditioned synthetic positive-control branch.

    04_controls
        Interval-composition and operational-observation robustness, Steps 14-16.

    05_plotting
        Read-only publication plotting from finalized analysis outputs.

PURPOSE OF THIS FILE
--------------------
This is the ONLY orchestrator file needed under ``code/``. It contains the
validated orchestration engines for:

    core
    pseudo
    validation
    positive controls
    controls
    figures

embedded inside this file. No separate ``orchestrator_core.py``,
``orchestrator_pseudo.py``, etc. are required at runtime.

The embedded engines are the previously validated standalone orchestrators,
packaged verbatim and executed in isolated Python namespaces. This preserves
their existing prompts, resume/rewrite/repair policies, path validation,
scientific contracts and provenance behavior while keeping the repository root
clean.

EXPECTED REPOSITORY LAYOUT
--------------------------

    code/
        orchestrator_clonodynamics.py
        01_core/
        02_pseudo_reference/
        03_validation/
        04_controls/
        05_plotting/

No scientific analysis code is duplicated here; only orchestration logic is
embedded. The actual analysis and plotting scripts remain in the five folders
above.

INTERACTIVE FLOW
----------------
Run:

    python3 code/orchestrator_clonodynamics.py

The master first explains ClonoDynamics and asks:

    [1] Analysis
    [2] Figures
    [3] Quit

ANALYSIS mode then explains the available blocks and their prerequisites:

    core
        Needs the genuine longitudinal repertoire dataset.

    pseudo
        Needs the 12-measure pseudo dataset and an output destination.

    validation
        Needs completed core and pseudo-reference results.

    positive_controls
        Needs the corrected longitudinal dataset plus finalized core and
        validation outputs.

    controls
        Needs finalized longitudinal/validation outputs through Step 13.

One or more blocks may be selected. Each block then runs the exact embedded
validated orchestration engine for that block.

FIGURES mode runs the embedded figure orchestrator, which interactively asks
which publication-figure blocks should be generated and requests only the paths
needed by the selected blocks.

SAFETY / REPRODUCIBILITY
------------------------
- Scientific analysis is performed only by the canonical scripts in 01-04.
- Figure generation is read-only with respect to analysis outputs.
- Existing child-orchestrator safeguards are preserved.
- Embedded engine source hashes are verified before execution.
- The same Python interpreter used to launch this master is used for every
  embedded block.
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import sys
import types
import zlib
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

MASTER_VERSION = "v2.2-semantic-positive-control-reuse-2026-09-23"

EMBEDDED_ENGINE_VERSIONS = {'core': 'v4-explicit-output-root-2026-09-23', 'pseudo': 'v5-upstream-repair-resume-2026-09-23', 'validation': 'v2-branch-aware-pseudo-dependencies-2026-09-23', 'positive_controls': 'v4-semantic-calibration-compatibility-2026-09-23', 'controls': 'v1-final-controls-2026-09-23', 'figures': '1.1.0-interactive-block-selection-2026-09-23'}

EMBEDDED_ENGINE_SHA256 = {'core': '4d1d28cd3293ed497a9bec3dcb10c22462d67eebf7cabceeb33a3128ff7b1116', 'pseudo': '9b8e88caa5b2454745e8a6e680c0c745b8c0be47c9e76e780c62ce20f29dc567', 'validation': '4abf7a5e4b44ec3f532a59974abfead9086f2340eb918489be0e214fe2550187', 'positive_controls': 'aca18f5570f1e6d63c82c8ca51beeb1b883cc18f1f4f448575cfcf0ea8bd3c3b', 'controls': '3c7048caa1c4c6e1e432818e9e9439df6489a2818ed420c1bef37d9c15866e04', 'figures': '510b2f3ced3a47ff4dd494ad371fd250822a6f2899d51ce844289e2e317bc6bf'}

EMBEDDED_ENGINE_PAYLOADS = {
    'core': (
        'c-rNC>vr2lli+_oMSsZrAQBR~*z$bCn8_-ZoalVEq$A0h&B_`WB0&i;2rvLBiB@#ZJj6W4yv4rBKFL;9ccU+8kdz~xojI#>Op`!gs;jH(*45SD{B2`i78{pI'
        'wh?DHY__Pb^Xx_E8`l4}&qnz;$*u-$UQPOc$0r@n^E!Dkx{k}Lh^o8@M|lwkvqfj<`0I2IvMMg3QI*`pJvNDwv_FZ;is>z|3H*DL=2^bC$f9X7D%oij&q}t{'
        'f5q~QRo8Ln=l#PE2Z#IY-O=IC2WKDl4i0zUvAx~1-P8Ruc69jex1iHG{%~@9bh_W^^S{o?oa(Y!k&ov90W^iMvRN{V(<F=GWmLhV>?SU%lFg!|h{p`zxJs%y'
        '3_D8M*_)HgxQaSOJcBYh6kt_!nZ_kcvMQ(6M^Tn%$p{Kh=V=8cW<^}aSrt)-Fv6@%=<7JDqB5>nl#M$H^f!*nQIX6LGaKc6Hp>fm2mMxMz|O9d(!>#fSw%@!'
        'vMgt_GM<m~z7A>U5;_`PW8*A~VL;Gx96>KR08XIbXg1!mQ5uycz#o<DGRf2Yia_TukK&E!lKuy5?*Gr@{WoX(d+cN%uzGNIbn+W}xBJ`Chck&_$6x0RFk--+'
        '5US#A45$T+r7&7H1}u%LbiwBExL}_?0hPrYpFRx+oeqOP`bppL&q-|Z<>nSWQl;^4tFPuO99>69e{qp~Aw1}{@o&5Dn743*<XA0caah8WI83rhTo7ig>M!~!'
        'yw5684)7o5fWzX3w=;#Iq~SP;uClzWfVA!UFP+`PXnzM3<;5bz>6%`q3%lZuG5?_=wALE0`dT4DaVyc2<w+U;Bjkg`mz|#$K)RC%pte^Rr)QuoTqf6muL&+U'
        '7+0R&!1Pawd@3hf%-g3=L1!1}VM5rg1ZBXmaapm;m?4;H0+nx*>KZEd`&hjX<*;h!;P7n!WcSV4!LR%5;PCjv*{Op#J6i!ej-a|&u~%_6$K`LVbZv!;g&K6W'
        '1GWd02k&us8~MDN&2<=&ZxIU+!vvI+LFYxlj<iLmQ3Di)QrCGpW>Gr3j#yVBax>W6@nJl;d`8zf5CUuT1}kY^0w(qjPGDgI$76@NgtZUrC_x;n;>c!zSA6<J'
        'U4`F%`UJ?lG@2VBi-^Yo+a+k=NeQo@l_bMyy@@m6DKYVs84}gF0~kcs%(ht<=YpF4^vS0p^x`Fd0aMSHOI@VLE@3nY3@DG0i2<)mf$BPlKQYu9%n`EjgX8^o'
        '$m@ST***Tnnopw5BD!szyh#Gwxu)Ph>u=xn`#*3Hw+3QPV4zw>dm$`sV|zg2qF>J9Q8GzJjIWPAF#%RH3u7gjZ8h6`F<^>Mw1cxwQ!fW?!iaN~L_^W)cAM?*'
        '3~&(vY$J0X6n|OLzWD#&R+jj4WeI?Y0P1hE8K4aSf)#mq{J(&EfXrasd?pqlQ1dLhA}siSY2>d4+8}GO-I%Nnr!=)LUJz5%MZPxITz1DG_r)3lVPDRtuuf3*'
        'Q3{3_Dkn11Luyefs_v+&qA7p?mR`-WdmV`iQ3Top;XKPhNNYGG{2{*q<&SM`Zf>&opnt^dO;n`7bYKCN$bLF!XB%g4H{OC`@fKGe)Fy#Hg_dA-McLvO#Q^e?'
        'TBCy%4mt#|YIb(AzppRv&cTG(X&I+5CS+!c`9iZ&JJ>)zu~Fq4g7<}n929JLqtoGBg-fH`BB@}3D$OH$Y;(23u^CzlNj6EYg5OJRiP~iVtFv$#Wyu5{1*7t&'
        'E(!WqY3eRbZ3h3ag{=TeH{Jw6;4J-$N_W;$iBx}&REAW=$Erc9=GqF7s<K-3QDxUsdNxq<0^z37B1b7DbXLwhfLccR6sRP@1xI=aQt3Fq%}P|Vr|c$5lQAb#'
        'l)9*@7kO$Zu4RwiUMHZSPa`0{vET~-(%*Vji|ho~JJ(RkXj&^2og&Xc7wMerpMH41&o<aw`1T7sK6-cX<~O%0-dIQtT_Q>W@bN-f6}k%3qtOeajp5XFa#YYt'
        'v0|{;RLM3m`dV)yIo{(MP7id#1QCfc65hdbzAA|;q8i7pqY^OGLtxeFDL*GO%>cYbbqM<WV@~FozN*zIrn5Az;_o?p6I8P5=78|v71@Y{6k>W?LGa4PN+KZP'
        'xWJW3fXXLMdxY`jzabYl9|}NQ%r%xL;EFDR_UdNhJ5Yv@nvmiE$;o7a%vOvq5F%s;+Tt^WCL@p>$3F4?+r0R!<}_d@G%ut%7cr=cHC{k(tJgpu^%j^Suz0^`'
        '5yurnXlL9m<k}mo;-UhLQR;g;sGI(FeDDUW?hmIH*H*+Ld8+Usa6`^%zjY1)l;kUFvX3<e`+cOiK9bvVE^gG}2u2HHJ=Skj=4v}?>&D{+=u;hmIjqyiH4Eqi'
        'M#85Lg&~j=yNpJkZ=+&d_BECQ>jz%u6KX+$hXxKhTRLPnh2dmg&4GD`Bwo>1kO5PUxU$ue&x$LU+A@~EQPzmBxYp9-QvM`?CBEmS{5AtF4pz-n{x0g4pmKbt'
        'e$1;Rt$)ogL8%zUW!+4<C^^)KSm|+)AeQ4~1Pi5pRhNvb1jttulfPj_!RMDei+MSmG8l#|hISo)0rCmgk9fXEeY{I>NgfF;u*Xh;i=f=<v9md>j1H(GCvSe)'
        'KRr9yJv%xHf89ShJvcgKL+0JQ?8kRtd5)4wP6^Iwe|vNLRe$qe`r9vHV!mNV8LhG@aDMb$LI1@3fu_kgTw!U#D*uL^a+w6mR)#bIPr2$eW>;U5Szmir#dno('
        'KLwpPNAKSs9frrdXJ7#%$X(8xy?^Wec>Efs1jXkMAIIPMAwEDP(yu=L0uW$m(GTZa+ZVnkDt3Pyd<@{DDsTR*?y5h5c6@qv^~L}3PkIHk1qiJEoc;0V1WZTJ'
        '?f#64BER+hAOFeg_`tKy_K(BUHzx<jXQu$uJv9vmm$%9#<MpU0*j3)+9WtUa7*O6<xMM_Zu$jEaIAcULFonFeTrr|N7!cmNju@!&0G#sZ;B-HHf3)}E-9D10'
        '?UC^W+V>7l!o%J7`;S2U*JJ%Jp!{p_{U~HV7SkNfuoQ<+`3LR&6qTU%AgcrbKtw!s51kHJbu7zo!=%i+J_92aq~$;;X&^1KnvDifsw+8ZFt0{_5S8+G*AJip'
        '<oi)Yh0v`r^vW2R*|_w4(NTFFZSTAa(ZkY3HasBiYj%4L#3v+8z{_{J!A*KId3f?!oOQk8((_R!UIV?P@t_u_$TmmU^XxMkjR{OrH;ty3<7gmC2B`9Nw>Gz5'
        'vTqswhn^+R^9P0%bREoRIP@+x<?9%^4zA<7adHKUMOTe1PSfyrFoL_t&7q?~Ei2ApY!~Jrt74(L6U~x*45$gL(Y^C&D0esv(dGq*BfcBOvx@CjRgqlIt9Tzc'
        'Xfa^y8+J^}4}1M0_<L<~Anh#J1k4cB9qe8`=J7sx_kYZD7&Gs_*}oiZoRz3Z$Ze2NxNW<rbHPF_w0T7_io@yAqYFt|+@U+LUnAprpd!<xENR~eCGfd?d43@j'
        'bN~!()E){=?j6wNMgX%uExX!`l7=^A5<akju%3ns4$t|-!&dL{<H5rYv0_~;nULF&XfP=`8Czr{zz(7r(CC<m+~*(@(+!(|=*-zb&N=7Poc$3ns)>#Uj(X|>'
        '$^u9g7GOam64jc)0~%*qqNw<o3G$eDpr*2Wtoh&({%I(Hc~M<?#fk{}U~d7C!@PQ%&$BTRXxBiLcOvaEbKg$e3?N4=nQ5X6?*CT~`~g5`H2R$c-Fgbc_?YeA'
        '&A2VAZ1i>;VV2ItYq>wfJj+Ys3DXw=Q^qE3pIP$td+w=5Kh#Z}F8G$jLK{e++V(p-tC6(uEW$-IyI>yCW_%|30&_Z@RUzu@M0#A%;F{(%cY6&p&CZIs#>LS?'
        'KjdI_Vfy6mi&GyK#urFyKY)o$W?kQ2v(TKmP-g|$j$Lop@HudEi_>M``IcKPtCDOU+dKs!>q7H>0O~%>KA6pg*1y5&wQ96p#q=V-1wJ@L4EUOPPh7sSj$GrI'
        'c!VGv=>!W7<VsyJ28*KsDcTH8w!degiS}QV?_k{y1zi&lGl+xiV!4|&V3RbDDnJu$YazAJgZ0S8fM5c|v+icv`MvD_>fvIbE%^1xN{o-7m7`fg3^!3aj{$IA'
        '9;wKl{$H^~>%4n>-)%ZCE8OeM=2P5K3D%tUE?}3v7C^jaTfwHks*@5hS56iT>v5!Ey)|c@Ix?r7TgF6q_1sHqyO#QgpkZ1p!}GjjYN6gjjqlKAMiX1}Qc`@z'
        'ZIY(aPpR!9pl76woXYuRlHB3Mc#0|>E<D7X+;J6)0D)t13KGxroaBfUsBHoRz$wNokeaIbOz(N050BzvCO8PnpyA7xljz<MZew@OdDQrq=vI*{V8F@fR`5p`'
        'J2b{kZ{j=90zsn_;&m;Kb)@!$8L!Vbke`TmJ-s|WiEhas$zyn)K^Fgy+GExQP3lXW35i3GuzSJecd>ttRIt9IYIF@c4r!L<06qKSjc^ee&*oVIMvr_IqA3`n'
        'zt|`wy~@ouN!IqxcpfHWctsvKZaf4*aKY8tY|Vnx=q{9-_#Dg<s8jWaqYcY=7Li^!K+{b<|D0s5b~$GY(aT0fx@crjRpMSXUqgG4qc_cwPhhbCzSau^$R2(M'
        'ps<%{^>(`+?)iB=M*jy_Up(J$*pD6-X{%m|>J_g}U+*gUe_M-HblKXX@Z?bwvp<~GWs$Zfga?qCATOY^6<JExi;>d_bkkx9fx(1Alk3)cCO1QSxMETL4$=Is'
        'XTwSEq9n`-m;p#h!c!ru3+?=Kyr_*=P6R+UwCP137Y%FnazW(itOBJ02(7SAVyJQr@w4kI%U4=nWnO?iqsEJNZ$TVUL*p7oYlG9#n&(^@M;_2`GG+s$Rf9kC'
        'vQFt!x2wrj;tM1jt1eXyu8Mp<>w4l9kbqK`(3|fliqG|O!Wr!<DqiWOVe-(~E&AiT13lS>-_gg0MFG!zx=lTgY+R+BbsNId+nZC1RZQA|7AJ(Qf}P5ZQJhb&'
        'QF$WKF}xRF9yp-Gd;EN;;T_(KFAoi4n76~&J%N7?pBE3UWw-z|-f5^%E!Zv*LD&IMsL%+yCT9H}r(QpZ#^Y{{<+^1cb7)iJnQkh>SLO_;RYEP}H|5eminnGB'
        '^g|R3=+!u%Q3rf2tyoB;f~bXaaX#4Eyy&b)2Hs)b2rT34B_0c%Rk|c}b5`sFiN)8LS2<sjCp|O^?A5Bx#|#SB<8l;z3sSb?piyCuyGxjK3re*4wueH>S!7Qf'
        'eVR^}6oOGD1(7O-L)5e)ymN@z`>h_^egO0W%^;azQQnhpiSxhI4U)sMu5_}I9S$Bb$u*F$sIMu~sxK%Vo{=KMZZV39Z1X?{Knm5nwAR0;s`BtJ^?SRvb|}6>'
        'EVK>cleBlzbR+>bq%u~Ij7uSsKjd9SA*9P6bgt@aSa^Ufv^ZMU*+AXG`GqgNkHio-T&Qi)TGWKj9p5?t@yuWah#!VupB{-xZgfuOWocRuX)N1eKMdJ6y;2Gw'
        '0s$~!&jb$<3x|v>y3N`}Dhq?PyEb9GkyO!63p}PZ%k@ZlDI>3=24z&YIhZYIN{%Xx=x~ZC{;pI(e6+U@ttaVN52|bu$T!j?CR*FZ2w{dMZNHC9p${rB$Go<6'
        '%-$z=@tAK4Y``;mDH<#Dc`@R1gsos)h#k!RJLh=?sdZ?gXY6XuH>hW57(=&#L_8_A*qCDX1ShGc=JU1ZMT5RSj?-kCRARpz8r7bl*rs7?`|W)0tZ9+L9I+40'
        'L?Al8(2#j`*<!&i1q$Yw;WgeiW<|F?&;{S~+7Q9L4j)0m)^%e`u~mb&ZMh*X4w907Pb(zNI(jYD50fc`k}C*a`(QNCY6DV>qKZ+&H2SERd}Nff%|j}w9%^M9'
        'R^YHgrMo-v?#qYpz7l&F7_X&sH?0MSCAqu}Cme1wbIA2n&Qr9~rqMZ+Kdk~Z`2rRp6V-EHH0M8$i-LSratCFj##@G=JI}LuK%uXwyp&zxDGSlbXKyi*G&8ms'
        'X))Y)^pMs8)0+VJYafsyW%$mb<|Y%{i*<a2)Gesig8)=|(u7ZQ4RF=grT`iRV#|1a_o%i8ngpdL&bKaLV4gE(g9ogTD`8b7+(r|T7bOUj^(xtM5h`REgT=J^'
        '7wYSFa|GMeu80Z2#V{>HsnxmQu*$<?{bRJGAb?9@DRC&nPB^DhE)8cj9-LxDnzRjYX@4N;kKM?4(L}}zYD0z5$YIJzdmL*U;A_-fXNq1np?XOzFO9uP;*!E4'
        'TcP7^IAV3yfoZ1+(+(#CcWh7D3dW8E!;S*Nxd!$cOvjI?vj*s|nxMa;_BOZlU?dq*TFyaHDHcwG=JmF=Lfl`#2#pYfTM}f)AvwmD)!IhuaXtc$0+_;6AHW~t'
        'Y_wQe(-9H1V%U;LitQdSr={fJ>IrAg$Q0?P;GhFHujFOE?Fo)ZTvwxMI(af+mwBG{I>d4D>!9I16P=u6-xqFwncQ}8?(w~iK>cd$RMegVap(;mA>-)co=#Yz'
        'Uif;*HXqlZW0QmX63~^QycnkE!ggHP`$n`^OlE(>BjkifRO>p=ih$%7#+Ck8V}QIZP^5(wIxDtJtg|d_=wu7-JGx5q%dYorP*)7H*%#dbNb1oVqZSWvk3^8~'
        'D2<~c1a^0m<nwZspjQr2!)e<h-)r-A(7O9&!puXcgJN10ajZug;jUsxrw@j<zT@FzkAwmgsvYxn4M$gIH+4zP9W6t!6jSq-oN>~agSYaSlPG-Ah~%=yWkho='
        'M;pDGXK9js?mmx_{Due1(Lp+5uy~&Dd6Z+a??VxtS9W#AXvX64XZFrWeHh22UjCf6_eCT|&#+2q6Q!qbm}pP3@OLnHG-l2pz^ciyzFxOp+1lYs$+*n|C3iw<'
        't!nRL+>5DSh@&S&OFMdNbfyis;nlprmNfZbNg3Rp6BNe}laN+kbZIBMT;L*!M|1Lp_-fDEU#-E{U~CzVFW!l^Ec$b8kzmhSxS$^A2nx8bSEJoyEN4}(r`;@_'
        '^S}$MaM@^8PP^AMRE=_i7PMmGiCDyH{M%17{%yzj`)iM19<_1dVB3L%{$Iy|C`C(4D}yR@m&X9nE8pxqo2RMHHsw~CQ3wq-{YjeNng!^r)HEI^AbI;dR9i35'
        'RF7#m^sqidl9>Fv%%WL+E<pm@$K`Gz6#X5X57<Ma(tW`ktXd}L@C`EhhV4_p&p0VZAgJTMI(nSn#6=p-a8I({VBrUR@F3p2pga-c<no4S%m~?3M{tmhYe6!#'
        '`(2At)lOj-^DL-o{KXT|_zQ~?DX7RnuWMPh4F(`iiZ$A^d!^{A{^bZL!Jx9*wXNxCm)Ah@5&8lai^#K<93fIPCnm})(n4xq%1WAH8Y6l6L`Z&V4HppECqzCO'
        'qFqR~dlmEXG;32s*duwEm3m`+KE)3+C)rhVqIRBWqIN8l>D~qhejo3oXv|2X)ivD=Unf5Rd;PI7(Oc#uPNO@kz|LPJ3t6GaZ~GKA(ML1Y0AUmVtXgup_fh|H'
        'uzOQMvnWN?rsFk@9(Nl3;ZbELr{dEQ5Rrvg%Cn_BIw<<naZAnRjvrI#pk?7yS`9Dr4k(D2ugWx)&zlh<Iv&iSt;KW5cJ_O&b~Al`xbEE2x;0f|v)^j<!rnAY'
        'CYGLUAwcZ+di6<gPt$&1JrOT>)#L?o6+X?PZDHGXzc$;G6k8;ymbh5=+H13&cw0Sj(p{(%ub48^0_vE_KGc+n0-kvk#n^KoBgXLIlUi>$vPlUtoX8}x4Hqs^'
        'i{Tk1V3gDFIjr-H2FQc4H$*N9EnC&In(Ck;3TKF6tUt;ycUksLW1XD11A8}En8N(LUgP-Pe+>(Mg`rO1=A>>3E1f3+!2p)e@nX<8Wy;Sn!mFNg!f^!1eY}Sj'
        'O80c{l*jSi4XNcY^vSK+a~BRlwWZ}fV6FrKS}jV0g-O&x(D?@HEUs2W?gZc-4ixdkAt*c91bpOKYM5g@&!bsTFbWXb1Mck9)p-?0OXTPc+6`!(JViOf-Jung'
        '6f26!pL!83%cODmo2#ILvW3>UMm1P*Nk)tj5jG3u85dilrSqQKAzIjA!GLWuaGKY=6q(qAhG|{v%N*Dcs0}l8;D&!r;_)R^LNBrow+1VP--kT;6$WQum2xY5'
        '<(MCyDS`+Pi7%Ayh5<9PIapI_#uY$mzQbcaZXp9Iid87<Mwfp27HW1tCAkVg7F{8Io^SQo&P5Z)#6>(>=+YU`Qb@8IK`%}}fY!N_4SAe%be%7oZDe7CS3wJ='
        '3g_09TD3wRpF^SQwfeLKcG9B<%%MA5m1^E7hiE+1OLN^SRQ^hV&=sD{W($g&VJ-11D8crxy8ij5KG`w>43w2$1Np_^YMcj$oLhSXnaW9PNggjB$EcdX3QK6r'
        '_h$4J=^{vWM79R26Ey5y^qN(+msjaI(qy=rc|pxkgBvTo0<H#L{?!I{yn)NOptd*u$_~9k)rQmd*s%B58nxS2Mx!8)WK%9YeHD$95((7agm<ChRgKG?jg`8-'
        '*1&_gOAVrC9O7P;our~g|7c{U0O|p1rNq}mj}l+IQ&k9s1G@^a(=IO{(=5omXd%~84aG4}jqC{m9-*(rx6XI=u4EFYV>1TVR5g|<25XD%HKf;aX|+VwtSJfB'
        'V)m;zb1UaEJOtCE<N^MKSHdu9pZptzZxYgBWxxyuS8>&)CqDbz5a&&Ma-rpC5N)m<;`$*gqe(1mn}Kjy^IdBS$>sTyF-mv9>=X)g9E)QARmZ+vhC*F27nhd*'
        '4Tg!q!;zMaS@$uUp~uf)+3h7DkfH_HbSj@1;}Y)qc3|nWDOZxd-lBV)_4d_Pm$m+3iNRXBj8lv<rafeeG<NBR6ui;%1fQm3fgSB^@yCpgjSM~ZE&J!c`}!+0'
        'y9XL6=3y6=IphEQ`+pD_jdTOuGa7SsAQKZH$1-Hfz@03jOE93J<Qt7oJV17x4Z><kJ3?bu5;^`T5W%XXPq$r9&q=@!(yMx$9Iost??v=><4h4&p-HnNyzLf#'
        'YMe4CRi~l_`NFPy+Bj^|tlJv-Z8z3>G(v044E45KI8Y~pcugIDq6Ir#yfFtdi8mh;p6{4OcX}J{@NKglkG2wb7>FfCVo3vKtH0wQGCCEY)>aL7y&%mFl;Amb'
        'Do9jWG1u=HbOQnwbtcxmZ!m)go_*Q(Rv0c`t$5kTq}uSSMd_fhiErunf|Mq_mLsg5j?^g$GsPnkP!Z2md&|0GKd`OMM(PEL=zrZkIXpQ0SwBKnlb8VS;U4?a'
        '&;5%TWg^czOp-jYC)S-HR$5)(%Y$4K%&cnlf($0`Y@E62iM?h!0}BRIo{Pb7NuEXp<{nc$qW+m=eR)(vQFBtV*E`?a$5q4;Axeu;R``kAve+V!F!JIZek^I8'
        ')vhA`hH6}bgvTpBs^(E5n@e*ENq+D;s;jFX%r5fVQe^mQsH+z)Z57o<OJhZ4t*M-u>ReRKEd?(sgvU*LiXI~{LG8y}Nk^ZR8D4SKP=rMsmi`F9Eu&gfz39{t'
        'Xt4MZhttm(tjT|rl#ZdtlU0QIXDuI)ZDHIxEQ}uNk7zFsJ-4+%6H#FHvGOKBeE%yw5JpE&!|+gj^{wgw{dXNlV<U|R>bNkQsNvG28WgA3gtuMpriE{$?UG{_'
        'K^9)xT1%$>A0!!}=Nto^t60PgyYuWB!NGEp4Db2vSpOMkH@dCFn;_ByS{@j!JaW!z-RPg_P;E0qpl?csbLi)S-OJt{B=exTuVGj=|Ce*ryYbN$5dc)u*$P*P'
        'hnjXA98+_~WTh~W8bZ^?2GkA7e|u}c0W1#TpA{}%jY0b#>fIe2?>EZCMWL5DJ=;6_aAuZ5fjKlisB{8#KEaLSp|%sHC=xewF(0vLmII?F?J;PVZgb+7mT1LY'
        'iDj@PB7cm_%*0%N#c7<%__mc#jy=)*ZG?_9N!p|9JTK#rXb5#TBy~J`(_Lyb2pR|IN5Hel4AYGO@m!jK+Ad0+g%ZS)UYnww%xdlSruTu)iKYuW=9y3MyeKsS'
        '$mF4Q)3}y@8I+mMs)ZG}Yo_(V0B+;0I)rbgJLGLS5`|#@UdJ`G*1~4l#;8SM{@k#6FsswPm@qw$6W*E?*XbSj`TKS|05lCeMJ4ti$posC;am6T$YdIg`r8-6'
        '2IuufwtDUUGnwz68j*?VldMd}ct{nm3FDcPXd}mN{lx{(I$v|6Rk`|3>B<1qB4K2DGAkp)jSM1)NMzi5s3vaVhd?!c^B$bXB<k64wX=;i`P|}{wkEXMw_A9D'
        '75R|y3x51TAM>PYYmOJUwaDW%23rzX1xj7aKG4I(@q$*xb8KwCJ9sW12PT9MLF+8L@NuJMXfs5wk*ZDRkq-`0iqE1^EHmuDLsUZQ0=vZ(OlMXd)0?iAg;ZR1'
        'lV?Lacj!D@aoyk2M3F){ckRmX*t?3{-5^1R4Jh|`dV}(pIGw$fg?qYz#7BrpNf2dJsoMxrwp2%SUJMm$=fkIthKN!16qJbB(EEAS{bTA%hofo24{z%?a`5DP'
        'f!EVgZL9*@_U01En#XDt%qx>HJ|V1nV5=VKO&m@Kx2%SqN<sM56Jn94j7-+f5A&;Av^z<JDJ8BI3TzwAg!wVoKs1_KidB0JtcPEjxAv(Jwx=%bQ#p0bt1*{n'
        '*li4YPbQOmJ$AQk-H!Q#%<~H~6=&i$ODBHTLQo?Jp3qZAxaio~v=@qqOvg~h7*j;fv+8^Hdbo+s_4>xGX*4=I1&u2KMHDZinpZhb7EHeG$E}vu7<|x66Hr%&'
        '=U~Q80fV@hE~6Axk}-ZSXSgwP(HxUj$~`!Q7g%R9lEW08bL>ly=%}O`%UXj`agXjDponx6#n*oVy*Q3CDbaT6b#p|f@xjM(j9C-gFOMtL-8Z!S>Pyn*fSCd*'
        '5fow#$j43EviQz+V#7#Y=alSdBWPUIMrm2JduYWHQT{?aX~7+iUp)wOi*3mWwqL9qfDRfq{5YU!Be@QOd=|<f#|tsI`hGzg#4G)NKUn%Vz^|{=-M#{*N3IB='
        '4cp`J=%&uF0`MHG6!-SHRgqXXff18zy$+09bdN*0RQVqRp}46k%9ym8H8(Bh9Galg&GTyOmXx|^?Qt-Ui$x#WeH@U#N>h$!SLf%U(Ke92tFpenOKejLu~#N{'
        '6b1Q6L@r$Px^$lLqO=KVkGgVD@0s$wI8<^CgmsK;xCE1BB0P_U^B#ESGnqKk2zMyfHEwFxAP$&PBHE+pr^@<&`A;-dUzUK`lLa{qtn2dRZBvcKDQ?mWo)#$F'
        '(|`d7_0DKB2dH{Abld^ZX{L?VE0~QN$)lI&%@=+mEt#V;g97(x9p=?;t$lLOmI`+57z0&yLUW|XarB`hxw+A`0*yP7*#?_<YKxiTiXY_dadv7j6Ju7#7NBmd'
        'nI4!{e;g$oH~evr!=vNuvqWViZ;Lju$R7{QD*@;R0LRE>VG$TyH|w`#xmVfhH~kg4T&1<UB;~zlji|B9YNh#b_0dq_>OB5GuubpjnaCTwFmr>|(Uc3_Vt5|r'
        'P#@jV{O$uSNt<{(G)u3!*iCvcqVIX5R=9pm|6<uGIG{5wsE5sFcwQcrF2e%6uIq%8l55l@<O1S*)^<;35>kPQ=9+6=Z*^3CM?d2MAg;Qt$yX6aI1}ypw0l;R'
        'xxVeIhd>!TZm6m|xz<qs(KDC#kY53-9%3KENQFYRhC>Z0Fxpz0?f~${(xvQt?p1Ai0dS87Nnq*rjXhq&*H~{*Xz<c+zfp{9$yz!7!0_8(upq%Sls;)Ytiap7'
        'qz!%jx(A`jZOZI17=8uAavxZFin*19H<Vw#X4(Z655@m_nn?5=s7Y)bH8778(2k~7&~eP*{{e+H4@f!dn5P4ba{wBAuJZxDKHH%?tIpCiKAi@|2|KQlN=xgb'
        'Eu+>tH-_<?x;(EfMz^L5xo>tEk0OjbT>>r&AXn;h<7P2)qb^{$#i?y;M$j{6tQ8x5q!ZaG)3u9zT%p@4vRp%o-xDZuiM2Kr&C0aDIx&^sY3wGi8_CdJo<+)E'
        'ryZNx_)a<HeOJ>=f$-z`bXE%d-~|bws#L@6o)lwIIZBct9gNaQzYz>!61lLpht7maD_GLANHQ5-J6wJeOGL;I+fH}=|I6!)<LAR}&4jqO>#s3Ixnn+Dp3~L0'
        'VxaeHx5+GxvVoAWGNCMmMqzkpcdN(#*&j4g6{GUmk{(eXix~nV;Obf1($3U*kC;VYLm$_Rz0-dj95a=c*>kae5n~|6!y_|KJ{&fKGe2hJ80QylV&p6K;ta%e'
        'hfayCjlN}tjf({;EDfdSV|KoG@*4wCF4(WN2e5`T0Ck#4^_WlAZ6uElEYp>QfTO+_`tgh?D{pC+yS~5>x2%^%zi0XFkL>)njcjGiglHQ(u0tniYlhr^%;O&k'
        'eg5zNv_lr!)xhTPzG3^sb7h<>WkagaOiimqgfgM<<a;<JS~~Lf^u`9A47NimF9Z0J6u|38`D`IQLD|i@cl_JgFGq(T4uAUa_U-=3{vKU>>TP-A+FPMu^FxV3'
        'XMuVYMkwk6i~Pl-ctX2XdUH|Vb8cssX!s2r-5Ix?nYDpTb39vCU1Qqz?^<h1(?Nv8uB&-9Gb5VKERCdfT%6m=jkr1uHH=gnwAwmapj#W<re;8nEp?+VQ1g}7'
        '`ESDbR<SorgQbhS?%n4ny@6u@da&|9pxt@?l7g4<So-H=HjBqKhTNw$>%eAbDgkG#^YW5f`162mIttX0(EMz%{OM{-21pN<uyjW?8ALsBWS`~nK+-*4$GVB_'
        'D6x<Ky>yjtkQ%lD-*9b~8-m(Bce5g|gw*B74Cp1v%8KvN;I_)x@Zt#lHfT0yjq1U?Gg@E7G7k%U_nSe~DvaJ>AU(4U-+3P*2tkzL7t9d7^&)IRUF+SRhJxiv'
        'cgg2-ZQDR78+6)Ws3@{2%9}^cjAKf`w#48l9P*w|u|u~|#}rKIH{xEMM!dPkwlsR}*&)E(RFGeq0n|+c!x=!Wh1QFuYGB3C&*r5#!ACPW_PEig#xi8`n9L+f'
        '@Ia3`Yr|I-R`)ul&1qD(?{yQ$aQH6kYWD4B;jL`^X;`_5MZ+JVgFkRUDpAVOCHHZl`_bp-+EzGh;ZDc(QB=BCrks!P=VYFy3p!^*0n`OQ_f%(^)G|!AaHo2~'
        'steL#WnpME#hI#ePcD9>5xqllwM(t#5=%)d{-ya=Nhk5cDq_kd`K#Ay7fI@eRcbxC^c6BmlSks4mPw$qK39JfN7ompAke=Vi{7Q?I~-mh1>Wu+yfdT#$}&_f'
        '_&IT52L<djj(PmLI+C({aS3d+T|q2fJ1q{4-<aIzZ#DDU?lhSF`O!bjfz;>_u@?0jl8c3|Z@K}tF-f@f^W?{v6$tu54_b~5uQ(;lla<loP|bJNJQ&PL1xgu2'
        '^swagw6gR%fbcM+gbX2RyJ3jY5Md}tiaTdd7l7gEKCrMZ#Y%wb{eJ<$H9}q'
    ),
    'pseudo': (
        'c-q~4|8nCt(%}Dk3Wm9=NQor>n8{>TS?zbuSd+LiWBY7R@}(RXhZZTDI})iSDSNd3T(xztaQo*y!o9`4$v(+-HyR)Sf|Bg<<g4#g$yg$R2GD3Ux*Oe%FTdDY'
        '6#3R=l5WN6jW=JG*ID|k^`+PQy626vC`qpd-lClJUgDFM@B6JRA795snTKVT2lFCcL|K2nYz>{CR_oO?OSAoD8qSh&;pvsV3H%Rj#9Q0D3@FIsNu0;&IPSMv'
        'C-FQh5-7Ozrr|PMl&zlqX`NliCd{PplC+G|D2}``^QOs7>_tf)kIU)OTfpPod;gvqegA&YYI*RXU-a!?|NP(o`sY93Kd)hsW})q!V4UTu?*BAi?)GMjX$dXQ'
        '^SFr9GAxrU4aSI#aygHK0-nS{lFIQa=q8xwPB751zrv6JH9=Su@$7QCYy$9Y51PD6$_312IK`SNAQgc$)%+*dOtYkje-5V5&$O%x$=8la6PU1&wB5Uk(>TY`'
        'a``AOlB-lgG>dbQosT=a!Q~>EM!{kZ^ALx#z=a6$3aHc#?Ykbs15&#U^C&3e@pYPv0TJn9I)&+*&BJlIX;+Fep1;&lcp0W2%*qlwFICIa#kgGX>}P326`{ky'
        'O*l;=@q8X8xdhcZ0P@YDQ%Jc<I$xAUo$dwU&+{w--~fSWURJuz@(+_~cI$=32Uu|_5Izi_iG8Q<mDh4f;L?IVm0^NCq8B^6JZ`VZ7Wp{#0O4_7W>DTM!^>%0'
        'd{1z884fHGpeyP_1rrxGz7>#qiar2X0h7oeo(T}y9Ek$}gi|lOEaLn|fGe+|^mR6kyl^_d4!w31PeQ2eZTGjo?Err`u7(eTfc2C*Q3c4t)_7D<a*|vvM2m$t'
        '@OHPix2b@jN)IMC%4S{>0|SB205J0G`L|Regf1=%;4nrRDxo*Qv<6J;iXlK>lX*N%(pdGc)dEbP9=$txHS*q$_D<iOjNXn8&rTI7{fGK?(>mY2*lZp3{=(B1'
        '%;K<E<natvO3^Qi8@&iuH&L-sq2M-&V!Pz_cAzg6QKSFYB1_$a*|~U<RSqq`cmXUu`P(A)PGM#CfOKIYd6T3py7+Gi(~|+&=~ZyLdH#~W!}}M!PRnE#dmn*o'
        'L(@}W5M^wUiVnnh(<|ojIGH44*@^J#Dvz%KLJ!xTUJna+&(*hYO_;=ryjLITvxU6oyI>r`KE&Una;HBJa~N)^SG<nHo23VfsAqIbu!LD0CD0Qf*?5-Sz<9$6'
        'JoPHx<CvmvdMiHPR294xA8o={6$5J6y%5Ow&hh&EnRp4(N5ynK<dg(nE;Pl|J!I<cZmj;?tSNW?oO}7w)A%_x*J_=Hlek=ZTVDJzL2h`}16c+5wl~j!l9sLB'
        '>Q8IL)jUzH@Iv6V*Fd6En6IP5k#~Bw_v6UhJM;FuvxB#zzW3`jaDtmyFj;Oa%*$i~1XQ$;`-UFy!>brrC#;-kDZ17*wuu>tJ+wve&TxJTWBmExSkl!L=n8oL'
        'byza2MOtPHAV(mgG@^YYNnz=Ml!e_$<BwcWKs&kcPDa0;9Gr~^S`<&?(f}v>4Y2Y>N-((3SFj(oTCaH%i9gO`n0n}o<Pj18tgbGMw7gD$2twWn0ev@NGR3V*'
        '4(u4aa=dqNvUM^#efM^R{mEn5sDet22$1^g!P!qo@6JS*abCg{#A4(=%~XQhv8)ffAi%>#SIk17$1nocm5?xrqo^jT5n3CkF#Z_M=F^yX8VMS~EKHLLY`DmP'
        'J6&jJQ9vW1U>#)wBW?!y0vH`KM+DR*q!%BOd7%h|ID%qA*4EGhaRU)M2i6#)fE|HQ0RN;YRP+3c&8TZl88^FyY}e#ZEfeou<^lqMd(n%2i*?xdPU71<DOql5'
        'S+oJYg-*coOR$TxF#iB#JNozI(c$UAFC#&8@!uA((ZH}$m?9yUh0_*wGM|YBei_3|#1$pHf3M9(eWL=M0ztRlzlY(YY#j$I@dP)|;3Xjv);17BmJ`G3^=KJp'
        'b0Lo4%Fknz$7t(<u7?H4;b~7_xF8he<ro;npoP$9B%?vMXe3FBb4Y7VvLk!l+$Nyp!OBykN@x|f1VPlRQ?25Av1EEm2%>@KH-JDf^;@bL1eVk6z)LZ$TRfwP'
        'xh*k+2#}_}hl@m&N4~ZcZT3muY)x{IwSr)>C_$hJ0xy}(0RhkoXi~g4x8$?@3U<gMmcPd!sqok9u(+Nkm+~ixC;UAt<hMB}J22xj`5UUlAFCgWGMTDhi%XEd'
        '$8n)57fbcEjA!#n0_&0b4m)i)o}#){sLmHrG6sgMUa68XNg`ZK{)SZn|NcEo86qrq9E>b>3_rySP%UXj#q)h?;|;O?qq!h(SZF}_D2;c$vjt3btJOL>dG*ui'
        '^z3Bs?C2!;Wpr|SaCGPmJ^$w0o}AR4ps}8yxZdve?u*{`pL)B`V1f>VS4W4h4}Lrij!%yE-@Q5$)k(76ADtfjco>|HM*HSlAy)4l?!Ec#^x)K8pw)VHbTSI|'
        '4^9sEfNevMSmpLx$ETxr`$xvBx;_eh`*CzQ65~?kRtVC5>)kQnVYK%)_~G5boBdH$?Gr>+zjb;xI{r?9+FUCYHD3ll>>d8Rq9;vysjMz*K2|C9c_L0m|Kr`k'
        '$!H&7wP~yIA8zIFyQ1bp3>&h}Mz#!D*fw9l^G-|B<LiSrqtj{_4N8*))u=G7ETEJnS*lr45`9Avg=QM4nkIVHUgNi5FB`r6;mt@EckDTm@<5#7LL1VaL3mu%'
        'NoB#!jx1X>@+dk3B0C%X`&j`0<ILRo0(Bnz{aDgp;O72)>ue9l%%!Rt|A94FJ)+X{ztNuwh+@}qPT#{k&u#%ktM#Qfe98~N`w@8xD1q1U6m+rT$y;c(qId$b'
        'z-^EeS-V53ewh!5vOr{5<f+;%`cSGZH;n$G9C!L*A%C|!eE@))Tv$>eR9}aF5##(6euo<>uEX7LUj(>`wQ;u@&?eLM#@7p&^a6C}0csP-r?+>XdtcKZvsc&f'
        'l6T#{0N3u&3lu<T3+rbfC$|0k((j<4BXpjD5@N|djO{^7Af&{3dm7Fzqj10_`zX4$s~&+&eChk0fguH7_ZM?Sj40R9@paw5jz2~T2sdTB8g~q<;ICPdwqYrN'
        '#2pRPhIyV{B6GShhY6&-fQ^0Vox`K{#}0M-BX*eqTu`wxUsfIC4v=`!r^W;n(Y^TNIG&f@URmbJ<)Vy7K(=|UvJV7JRRxs4EixFS-n;}U1_r{)Y+eRZstc60'
        '4v=~oPX{<48W_&AXzcbHX5Tx@7uuZU;Vlf6?2_$?|9g5zweJU>-{}K{WX^Mpp9JZ(#46^L%yXcew!f!Y43Qd8XlJ23>iZp|iqtmML#i;O&bFa}P9Ma5z%^(@'
        'xz+;%=~^`+K46pC^d@d&x5woW$63|nId}~iILyk|*&>asX{bgr@v$2;)(I9Icz4*q{kkzQeyB4^fiu$XxZy(=K+jX?+Oz?t6fA2}fo8UCKXB4v&i6H=?<6Cg'
        '7m2RP+WSb_&>DGK>QEDJ^oC}ER!0FSK?QNcxzVv9>VrmGVt3Ogo;9+AOEdVun@qD188F$yaT4jlX0vGY0N3%+hF19VJK6r#{RPad0<`%Q^Z2N@ax^QJEU_`@'
        'xG0aaX)n@WutekP^!KRavqe#YKJTTVzQB6#Z@4(yfL-s;+=)MXJN@m>x<(4XTrpXCSWh4g>uose)DTZQsM`S$B8X^ln@=hmb|fKSel7k?v#UVcvPI3?NdFn9'
        'Hv`fRanI@^>jw}-Lok?u0R?3$6rg7xa1Dqb1u@I6Aon2qpz%)bBI-YDw?szW@9BSo&FmY;){9U&otD1Z&Y@oy-kt0nE~ZMRP<<A1#VrVZfF(M`9{O~O4Noo?'
        '6JSvIXosW@s8r1;tMt(IO26+P6F6)7R0q4}i`!_3gNDBy^9|s92>+W8i!#a<<xp?z;CNIk6X&^J=JaeIgl)4_8GjUenE4n+i=xqRr+QowB5F-&^<c1~Gm>*1'
        'qP0O6#qj(J=g=(fS`AndQ6CL`i|hf}eKPb=nT}2USXqL<4HFQ{I-f;Ye}zUjb}D#!LH)lgt4@28XMc}VULW*`+I~Y#50v{skY`y*A{al)W88A8C(=-%KVvIK'
        'rae+73ARO7Z|Dfg{P|KSaL%GH?ar{lcbvqqpOf84u#SGee<9hG*j57?kU(x*8ys>UVlkq2g@&%iiLlVc1#2}>4HlvT;oih+(i`>{s!{DNPd6axhqS7IX-E>x'
        'gn2(ga%}#2m=tkcG6d%IHYo};zfl*k-^S*TPvzz>ML*$L@|ZT%dA8Mxc&_*&vcqCvV;J>o)t0Xa$sb6tTb^k+tGI`6)0ke#X+6{G5*6wdD#Msn$|O#^ie?eB'
        '>(z@^rS4nu+|WzI88I`%t~D?O$l@FeT@VSvd~Q)40b$0I06QFsnk4>`kSgT1qA!#$kR&Z3pA_oOg(Lz7w0`0ZF6y7KAY7NZi5`H1OTfRUzESlX5B6JD94Ef#'
        '^`ua0jDc%+v^_%z1+-vc*wyN-dZrdCZGiYi+W{p@6PqLa(B4(djoyl%+IvNx;x$c;Mzj}NQJXmR8a2YzD-vbur3I<lMNG1^%UeWSt6VE@6sa-VVTqiwT~#7_'
        '>@?Fa12c>qQ|7S*;j?U}C!kmC?DjO{sTVGy<d5H3PsxhAhi!`m#b#Y&*Gy>%7G%8=+yH(eHXZvK!2oeF`PleYn8_}e1XY&NN|I=;?~+kcjI*0K54dyGo3`Fl'
        'wNc{1%{lSwf$Aw4nhZAOP-pd}TtV)m361udn%=7z6Dhv*NQ}Eb&2DkKMM#z5mFi#UR90kp8AomEVb}W*FNeZZ=LtoGAR>|_JZY68LXoJz$t~j2jf~plN^Gn{'
        'sM>I0nnedY-WPbhWHa4LjkML{W5XkkyD5}CcAvx3dU7${*#(N3&j}s6@0@ReWM4XAFamu>_3a)}eOGn=zVF*ym~pOQr0Q^0k(iB!hKkYf5*Y66R?pNbV5XyH'
        '%pfMwMUAWR;*gw5P{EyxfNEgwJ31g^9Wl776Wy9Uy3=Xw*ym9h_Jw<3b#MU!0q;TR&G325{b=A0%Hl>d!D}Kk@dv?R=&PyFOEJ1;ZFg$LX?%+s87gcpQ#^Jg'
        'MvZ=4-1M=X6ME1V<+|ugkj#+34E^6r-(WDJYJZZ%(}*Yx#S>I6Evt~d1}Fa6=`VPI0P+($xsi@iv!Kac4#wZEwbbc)=XlNz!pN<l7x;qYB1$TjOj8?3&xf{V'
        'dX&c==+QC;+U7$$v3<}~TxQu+p2DzD1{_+qUJAj{*tPXj$O~N|2&dC_BS$VnP>A6V9d>ws@wSa6Ae-;hP6}R!(}D$8!@8yzCg2OWHva?^AzM|XXAY(DsI;TW'
        'F^D1smVd4R+kC`ooQaI+4^GUAK{77yWM7jAPye*-LHmi=qo2fI{ChxrvAEQBAu)_5q-tA>hQ9g+Q|qIgtzbO4GPN-S1>ASW4j$LxPvxVk4a68Ejfzzs-4v@e'
        'S`E~gua(r3L_qk8>mZ%86cooPn#|&oWeaVoY13bU9w)d;gqm>Mzk+(sVf!6}lVB|z21+&1qRBt-Nk5FDHk2}LPSERqI7hKsyptHQ*EFuEfo5A4Z*i0p|J7YG'
        '*nSbQnv~Fb3)8E(y;J86pkwW@hQ}>I50D~cL8m#~pbf#(;<T+{uLC;_iWgtg3jir)DBhA{c4G<1f3lELRa-i@TP#f%kVy&~=mG_siCtW;t9i=@=mO7Ai5pgY'
        'ok?y5Z>?jlvqB{Em76T9_w27<>AD`yjB<bLOAE-8@;B*%83>2Fp%i>ozuN{;TRaEoxoc8AvmNoy6ySONXuZ~yyYV7NL%iJ$eLnl9OjXgnH3lF|1KOryanpoE'
        'Su705v#=arlRGEFC*imjRbooUz7Xr%QXbtp{b{cU!@x-<fF1Pwto#OzQ$br$7Q7W5Y(NYRqe7b|Yt~n|=5+dkRN{|i+m@R;+LYirPl1~PYo)18t<4N3xwq4y'
        'tpUqDg!VVWJ4#A;+BhLokmbV*zZPj}yIeZ(<{Akua;x{+@_+J<^J9d)L5BqoR33fYv701B-?J-dt*xJpX*Tm37gYLYcNNQg5Kr9jtj+iRA9(<Zt)*p#zU?u3'
        'a)$ZRgD>T>vh<RlbjlvRphV_sS);qQbV*mk6+w~Wd>6FTJf4v6IJy=K!0e4L;U{s?_f7=@`TQ-uL0*W?!sKa6@dwI<7+OFbTzX_HCt(@VkO*J<3Ht5X*zHRk'
        ')LuO7=#dm1>3Q2flbG9AjxO3MUks}0aN_sx^yZM?+)E1bE9F!w)Ag=NjT2R-ayHjQYk@lHg8=We6K?2PFwL$=0usdp2jZ41trED_mEW??91R1IXe)11E5)+t'
        '#~<VIf{Y8ApM!N0>Gu|Xw+2GdN~iwLX+{vBFxUm2I{<2f6C?BMl~WE1C>-^jT(S;Aa3^PNol9vIE7>d`qa<J-vs(!7+9#MKUJvCo*3W+xCYegsoCBNM@)%j{'
        'wo^pR$1VM!I$Dh5$Rr}hb49rt#csHY>;n>ARkb8;!k+2~pN|-lw~zs$V5QFvYB4He*EoUw3Dg^@`?I;FVqg&q|01?#;N1aK+m%Lm==aCje7W0(_g$~9hj%I*'
        '!sL+w(G&&nsc)2d2<nV3C=wA(ty<-#X>a?pn6OT1vycLEilGai$81HyhH~r;gsrn~Rxvk0@^%J1E3a&efYq}tPj{5o9{an!M$dC!F9)(G`KZV29gcH!?3=Qj'
        'IJ6U8avjjhiLw+nGDN9*SJUjW?Q<hqtI%Aq6MU`wl{?BHu8^jr^k%b^Wq)F)H<O$VpYkokr(BCB^IHyoM7BS%iODuU{n_ABHKT}OXp$P3t$L9m_aWvT%glBs'
        'hJ*2FA2sU-C>J~C&t2bX2dQiz8(VpldFRE8O-C>Ln2cjhRm1%WfSv0@xG;^B*wn+Z!_I?xn?tsbN{hvd?-ad>rOICamd3Jf#maz#8qw*I2|C_og$>tJv6xky'
        'V8?U|M@E6iE)t?xX|Mi`KBfVBtgqn8_X#AaHi6k<zXTC%?N`6GfNQ@~t%%2is9SJbx)P%fs6W*b4YKr!Kp}2UB9QPF)rE$5S`gHu>^4=AUSYIkCE^Fzte_KA'
        '@x==fUoY4>jUFHGb}^*h89C|^Itn4`aQw6iyEI}cHD2`oeIJw;+~s9XiYl<MkUoYYL{UaEil|~&u&o3PJpqL$;w|tQ2CYuQM5rRD*mwx&>e@mO`-+BB%BkVq'
        'hE=@m1Y5ZV7GfwbG#K_<YUabH0hsF4k9I3i&0M?dE{(<vxs9N#5n|w^j5Bi9ShaDIOvF1_n1{0hElLQV;&Hm#S@6@HVQ5HM!;|GUy?)PJ18A+xvYzWFo3qo6'
        '%s@8$SzLy7c@<_`MSh8Jz?)}cvUOgHb*%}ZM-G-`S(R19h{m4M5`&Qq%SwZ}gdCvRSADGkr5NYQysW#c)r`xq-hSc!Sxc5UPNW6(iwJI@$?TS6onh@?#uWZO'
        '$Mvg>B^2L^ldnheXL(d`6Ug2~ijOUdEu`S4SQ(aklKMpO6TBlV69Uj?xgZ;LSTU99vf^YIF(A>}6Gm(i(O+e|xXl*R2z$Ori*SNR;+ZF+aks3<?4r!(n4|{d'
        '-jFxKPKwr4bWfmBOPh{~IS#qj*_mP<!A7#>3ch2=3jtD|8{|4_mezU6)@${f%=eY&pX4=db97xt_quj_+Uj%)%dR8sHNJep9b{CC)yEU7FfU9J0c9IDs}LRi'
        'utD9?kB`RtmOY)Ob3bFT7xqA00dn7D<7#c6RZ^7+i*;0mH`BImqH~qeoIU!(C00a<SkVfZoB0Ju@OSr}mUR|)XS#dt27H1tFg!#|tAOAt=Z^_iu1rOm{YZos'
        'yIWR`)xwW;S3lpPs|HvW@oDvmtN6@U+w{xmAT7iIGa*5t7cUI=S9y{m1da=}Fh0b~E(f}hLzBQcdX3ErowP-G5X}Zp70nw|H^grEjhct-su`?c2l>!Ep~Y$F'
        'g_lKJ3`;fy%%2BJz3X&P`P_+nJKvgywMH$;YN6UusT%yS3|pdDHO?W2<=A#DC}-BJP%BTI9c%p#FmpfPooaX|oB1NY@5+<1;T<&ykDM!z1p#gG)p^(Y1xQ7q'
        '591HT3J-)L;DC%c@&|5_#imOw0Yf8&I{37YnbdSgb4u;yq~yw(yppv@^BQ<&2#i%IDzOpm^>`YG>0%yGS}Ya(u0y$RSRWepRxS~evRvyUDH2LrH;x5acLl0C'
        'pSn`S8u_oiMAQcK>1#ze*>jT{A61@RYYA{F8Ho=r)kZm{Bs3&Wf=EUOYp7=0o3Q_iMLX4%HYzs1r|5t03!i-YD*W=R{&woV`|pa_I20x29tOG07S<u(ZJ3g<'
        'o5%8WR3v4TX`pJ7f-mINsYsI7kPE7Q_F65mjm(rHfWQ>n^eh!MU%F~lEu?yjM>b~q47<#3ptD4fw*!!rXJNk7v{t?R^8m={CJtn3Kp`8%=}kP%=5bBcWv%Kd'
        '-1Rgpa{R!?sB~Hy%)}kd*wm~ALeO?2R4Zn=Xkkh!Sb89z%q^s9b@X&0(`>96FdD0v&9FB$8XcnY51WTirskKQi<+ZyEhS$jX$sr(e+nr-{i)S)ocDyALqmSW'
        '&vS;w$8_w#Agmzk7bN2{UAEhfqsJ=0K+_3p%@q9T_RSa(*|0N2#B?WU*73V5%IPNmRpoq3yGPPmV5%a^py{aDY#OkS_N%oFWWy*77zGZKI>;DR`jZNonhwL~'
        '>NVTkIwI-VfaJW^L9L(h8ZeHgv>~;d_tuU{2o@}{{cHrV+F;gz+z^_xdqPJ0IbI40jRa_YVu4ewIw_IUy4jh@=Ci$_FY6<w{*jOy=%86Hd(iB!5ll9n_U3(c'
        '-DGq26EWu^&FdzXXEpsCxhj(!t+hfcaeqMMA+;ekLaV$_zVGi9AG}*UdjcsN)B-wU%hK4&y^H>hc9LE;JBw2~JPSj^?ECk4*#7=K=fZ@a$fdI!v)##@m_^_7'
        'YMF73*ngP;70<kITPHIyO#wni>7P#6%#65XK_)}D2}@eSYdJ$<i;=c=^OUS>+1hzu5^g)nRognI6USdDZ-{d>)@rI;raViGw$~!*5HSmp$x=6V>kckmG?5~%'
        'Kn1wEiu1nrN_z`0If^s8P0@!hO2Vr&19BQe<03~|%|TK1X}FpVb%nvey?RyZ6Vm>AzW25F^0}Gm#7YYds4s;P^k4t@r^jixe6z|aZKfVsB?jKzSFeuV9={o#'
        'jlL2`Q!=qD5VQBy!J&K$O!)BV%==|@a`5_KwEvZsqhJN^YBa$=#N)&^Cx6=cf{=Fbh|=&OS}oLU>^H=|L(NB_4u`(VzAI+N?w0?{-kXDcugyp6MaLH*xf`8l'
        'nwA}$)V6<cND#NG+4Xgoe@(cj3jW)J)6;{)AN`J_^Hp*)c`l~rT>Ai>S)<gw^6D_E={8_|_1ox_$h4YRtFK>iU%y%npC%^8B==Qsb@*c99toUXA3O^I>(7qw'
        'w_d!cLAGr~E_;h=#10_7+xV#Ye(TO@w!L@w{a1UR)mQ)5y^}*g&_HE%7a83PcfxOBV3IB((v(!4D;Zo+<etYJP+Xz{{7ojGr_;=bB%j}ice&6n^2Cb`k~<-C'
        '#8+N2smHq|=Is7nrnCpXfyv>!w#~1$=b<GSW5;;}%#CTMo2kSx!Wr6{m<mrP7cx_}=%PxiuMo95t-ka#<+N+T_&Q6*$a(!^--$Gp3D4X<{8N9!6p^n)3o36M'
        'H%ZCkRL(dVCB8m)9Kkk0xhi8hjHCApas&j@!oMH<`+>_rS&pr44tG$Bn&5L@5*OZY*#b7*tRQ+MItDR>CxU#FIHi`tp}2_3cf-{qrle0m-`mq%w0K5m;gY<V'
        '=ro~(MG!|ApV6<w?YP{@1zKE3cINCWT+TW6G)!{w4V@3@TY_BGMs5S{{Bt3CUU*eJkijwBhG!5#r&F}F%e?cE^QimGb)y4S#nn(9d;<pa(wj`ftD<(w?ZOog'
        '+eu@L4O*509e4_SlgLor5X6W&;W=2XI`qjtP(ZCCSEDvi^x9({Tp{qAI<xi8K(0}2yW~5N^4Y+Fzj_u&t;eEUUjTtHFn{3e8?&nuW>jN5<C;z}nhO>Qt`cHJ'
        'c1Fa;%!uYR2rRFpT4#dBW7tw~fg&zW<vQ(p$Kr&9E~yX%-9;trL_DU(^p4g~!d8KIL6V6b@H(DAH&Cs$Oc}AYjkb!-VmG*bv;evcm6r0p!It{}t39z9ZI$zS'
        'GotF2nI^OwR=N$KI;=ITKx_88O@KcD?usrqn{<Q~2So5Zn+4PODjY8Z6sp?zPet)&)TE+)x-td($~#ci*)P4fz*#`3#{j>Sh11HVN?Z{H>rw1L!eA4o^_x;O'
        'fOMgx3M%)M{eaCeh#$=VgEuGjZJXB}XyIa)6%GXCb|#l>S6fG2Nh)1qcI7X>bMF|t?(<i-x80f2sN@q;7&6~h3)=l&=IyI?x%zuak*f>2-@AsVj<_9T`=CL2'
        '8r0VO0k9{X&vBl0B;-O`SK*ZJ3;1mLgG_#<9U9d#HLOm6T`XtQB>hki*=SjL#3;&Qlj<BPm(O6HVm%bd*kTF5A{CWCQ!ZGF@)YLB%E9t8&4@f|7i}7L>}Kz0'
        'GbilY5cLQ7T&KZbv&;Ge^!gGh?Ucz}@)|q9otJ6=OI29R*-~B&L#2#i9xkpK<_s`vLYG^nLE6ZEM_@E>E+ej)LXf35=9=7(Mr^n4rq=6WR?D6l;4D*ETJg|O'
        'SmN0c)x1eWGe8jDlA2#cJkQ&@51m;&?croOnh~?;EZyU&wQ8C3ZX<zoiXL&&o^Z8WQ6MVyw>wr*;d0j_mpfU}jw8`&s+Qc-`v#~0Ah8zLQ#Z(TJ?*KSf{Rv_'
        '0Qs@xQ%?qplg+|o4}DksMsl)*SrbMrc)niQkPM%g+Ez}K6<9Ap8f?IV$DSZP3YO^YJcpi57qhfi0Kv>cHOW>A>Q9-Kl_^(OEbC9eo!<XX>bgx(IFj|mWQp^S'
        'A*Yvr6gheRF8=cj|9Ou8eCt~b=`u<6hWy=avjoMg_WstJX1Dr_?f&+TyWn+lb#0dX)1S<evn1`MMn|{n1t|Qf*1H$=f*5!jPECOJAfhOnJkdhCF!BW4wOk^+'
        '{OU-&O@E!Kw)v&_;Yg}Ye^nfTA<G~sn3Bxuv|>*owzzZcG>Isi3ppZP@<2;)y-{dd^^|Yq7^o;P3bVKqmK-~YY2{b0()EgVq)YjI!)#|cx(AmsYdgtM0Rawz'
        'j!2p6E*#v%d4bo*4LsMakKI~Uc))&$9`WktG2Po=q8rU#sM+Z`w5p5cs&IgM001Y`h05YeOe2Mv-j+r*6t+qv*~s9vI|wS!;`al>Xe|0unY`A=L1^U3rS)7i'
        'E);4-&4pjb*=3SVvnxVvu^@6Qtgr?ySWUsII14-TB04nGn}AyE2D5R1hQ^XCYlPV(on)>@D%TX}oqD*v9rWQSXa{~6!5?&ky34y+-k$Dmp>k}HSN6Kg2`d$1'
        'mCvne#qJ%P+QQg8PtfN#Kysh)s`q6$M4%FaLKA?^`4pb*{QyFOZG)kQsD57ZqCqeutp(ftZ|Q##3WGkxNqS}1xAA!I*eM!%&W_yJAgaGv`rHOCOX@d8V)2hN'
        'AWGqPQSVP6Zd!nH=Q~|*x9dIYde6Jww-@$g=UIV{IyXU?FVZnHw&3q^p4CbF-iW0m6-$sOt5vnvawlP~)T~Ph88N7IU9Pcl0xt=n;^SD0fAtHXgR<Q@-@edm'
        '@U6ypH3om}%P$5$7n@HNjK~4b#->k3H$J$K9u~$8k(BAxPDstW;)P{;(vm{!g0wG=^1HG}!di!heg2m7w$#OZks=WT{cLHA4!2C)2C$G!+J)cQ_JUub-GUm@'
        'l+Q=pd1YtVzI{xFZNjUbVf$BlpfYS%s2XYsuC40gZY7Z{yX{OTYUbg!Yr@;M5LIM`QW8$E%PgsgWb3m1`_E*Odvowv&2Mx^crI=}lvjrL_wuU+?vKaxQl?bG'
        '?-*+5E@dRDB2tA#IJ6ym(6?-V2?$HZMg7JWT{;`c^pDUnEvEsW=vodYD*X=LjzG7p%Wu7xUKMvE&V^-}mwEtKSEmrzRVW4$uhJOS#=9@JuHvV7I)5z->(`A4'
        'RAC0@Amlo2L&*_I5=`|J<;G#90cf5nLu3-wQ?!x~eU-{lUb`To^Tf>yc+|7AD?{|^d1i%IR^NX}r#0#F6MC)RWLdi-6EC*%j>*M2j@X5GmXsKIy|$~By-@T?'
        's7^*_HXA1yIW5b>!$w)VplF@}`H88D;w|l}S8eL}Hg~eTJFai9g;DOyt1q3<$O%ecO?6V=3>!1rHLARu<YiP0&vlO%^YpE29#~y$hBE&RoqSRjR+2aIfD$ut'
        'T!WkhLO|d9RU?OPbxbL<s=~&^`$$v|=F}ig*?GKox#UnP+OqJ>a!dM!;2bCetKl5p_aB5k_Aco$LNr8*^VcD4;AlJtF-*s@SS#vcE9#()7%??^`0GY%q;A*}'
        '#EsT~YJ!4SM~ANuemo71PmcEAy*dMJcqLT9so8kr0Lac4r=!vSrXVFfH3Qu{+<Wue>A~rPz{*1c2fXdg=CtwPpjsQVcn0ecSHfjD#RD{PAHf_);#dz_sr3^W'
        '<;xzdHvs4$p_Mh4aj~kF%9e+fSwz6&{5HV28F)u?5lnrq?pvipERY9?wHnFrpau%Drhw6>FW94kKvps7r02uqeof)U+hSMc!LQS^y4&_ILvP2tV34n`W%~xU'
        'CvK#L*ZQyO#!uniDXtW@ue){iC0ujwt#ml8G6)+^@Xfc9`yAw9F!+@_d|{^jqw#6Y@~tz5E0<kp!z(5Y6zdmkHQeKPn3a1<b`hx4@29ejqAm;>H+I*CsaI|u'
        'r=AUX#f8{vbY^+WF+hL$fN>iE43*heBoeMJxfTRiqtc2(ZHjhLeRcWrQ&F+H;?S5P>cvB^Db=U10_(bam&R*Zc-4%ts;W(#yzq}L$&_z(#l@^0!@EExEElEv'
        '>Xf@?)<>Oa_1$fZGR?svrJV8TrffRk0_&pmjLJuKehs^z$dD#IT+vOQ;UODf=jzKCbSp}JO@RPmDTCi;Y`OyooB>YQRa7*W^uSOo+jPm$O#Vau>$uBOe(AdM'
        '6%X}Tu!>fYq=xiww|3Wz+yjQqyMum<8}7+f(VA(>G(u4Zd##|uIAz$KE5J?WXzh%CsM*idRiQ9BR$`&b2~FhT+@=@?<2DH2Z9nD;J5<4lG5Kaywi54r^zX-`'
        '!_$LbMi;)r$jV@M2`uZ-KrOdn9(l@m?xhQm>6;GCe3?>IM2sv7R5Kg2aNl7q<@(R7F1&<Qu^5*N!OorC;7nd-6ydNQMQGU`G;HG$jRw0qqA{uuU4e<>JkP%J'
        'Q{IX=FfK2&0u*&sxc-5v+M2QJ(w*JJBt3Us#nS3E<#o&rH7#yu*FSb%vRkiT-D7J*IGS^=cR;n+RIDC|v(Kw~&mFN~dc9^w;Xu8-RmF`21xi1uvZ3s$9Fs0M'
        'Iuw&@*l>){k6LWqLqqDUND~b_t;6pu?ojBF508kXYeIOd7O|z)4{u%j@i@WznNxYqFB-1UvWh$7vQqFxF3w8dMV$oH{hKhE5_{o9$Y0{UIT-6Gf~k{nqbcEg'
        'FT()oTNHiBCmPL~FpXjH0p<~U!0_U5e9e34Ed?~$f$J)&o%-5u>0*r2c~=RHEa{MOAm9=k=i{nt{#7FBn#n=Q6G>@ADcbha=Agxm6JXI!C;ZphryT{OChvBI'
        'rH+X+lJJzw$JOmigK#S~lF>ch<ypVw4#j#E%HEhNg+ftk;d!E*d7Qw#KZo+zWgqqSmNULbhi2hW6N77UQ}uW#oA}5`S1$$gS_@!tmcpjj<m6~(17=W@<GPbC'
        'aad$o&KF{la?mMfDxm0R+LCb`GH5}jZ+C!mBy{I##gwty{dkxT!pOVe045hYhAtJ^EN)r_b;Q~8ZUxM_qo*7-O;Vtznku9;-h7F9v#^`=dl_4F_@#G{dW#g#'
        'DS**NBC<fn4r<}I!j#HcHNsR2_GsB?wt_^uGK>WID5lxK9Dpv~!x6%Lo06M@ZR^xm!`g=AxQ-g1Y}%}7(s~`Md}CIZhC#a3Dfz)0$MLAL_JY)qp-@Js>3Wez'
        'X<vF}T71KHZ91*uN(xY^%cULPxS(N4X8o&-wc5*=lS-=KRZEloQmM^VqDSH+v+pH=A|YPUW0_bynI`4Z<K(X*h@OGj!LwbA@dlwlUKs}yHq~>bt9%UwXQ&n{'
        'ObRIEM``D0FDki#i9&8;7a|6YLZ)yrtW3GIf(ciox5@3Ot~F(*JJl>7=9^7v#VK8P>j@YXLKGF*D0%$O<Uq>0<o*v@h)&zaC#|qX6O6LQVJ*#uQ@Ntr(>8R^'
        ';U8$Q(Pn<Qh8!MfiFtZckkDv{Q_TamaXvcQRGj;YqIAw)AQgM$gvk#&mIY%ZnHst%r2|&%<no`HsWVMKsap5WeRXY|e)LqCC7eF<p8scc+orrb^e?u%dUb-h'
        '$&4%C?p4eG;;ChNmKcYgQ-qQl`0weRwkzJF4T6gI7xOt<uJ3(kHuqY*xsjg)?NG&aJe>xV$hexm&ITd$wd&F&OFisVXJ@(qSj$uMdwRZq@|$<^?(hQ7A5ie3'
        'i+e;emZ-~8Swvy}f{4)7qEo+#%DB^}CS!(aaZ9d2B1wZHVLcg3b>5Cd{cqm+Z(HfwyeABpFFL>TagOO+Rh(sVMiWXdQT*Tk<4(Ha?_JU@V|I=Ob8L^e9kT7Z'
        'EU~jmZFv_iW(|i$linmbG>i8n^6=r&xqtlI*-uA@?+$->_xkndWVBCw!`~5dWrioK<o(EZhd>o^WrV19ym$7KxJugwNN+BTLmMNA&sJ>BmdukFqaqM2j0Ou8'
        '@f?5DE6SYx)venN_nxr%(D2sB#cifc=ETn{d@RMS4lIB5p6bK!bi=(^-k+r*S~c4NyL!TiRQSOLlA~r-0zAny^LKVMI(mK%6i2gz)rI5f!25A@I1(WehkI{F'
        'wUBmlO692eIbK|~4A8b#k9}7@s{?CrUHks?%AeI4bz!@=4!mvarSf(iIO4y=l3jflY2B5cJZAB*;^%{7du-B+cA(wx&2Ee;SI_s#>=U87QRVCE>sbIgj$P^;'
        'o4oMkX7!bT{MCirpUmx~kuYSo0B1grrZg<w>gHkxlcnCV?segv5U9v@(xQ#MvGddiyhV5I{pIqbB%9%CVaF7BhA!D>W4Zp(;`;+6US`-$HzTY3iN4f&BA>5I'
        'hS}Zbq-2;33tm*|7X@WvHQ%~qnQCdi44<fFntc@gkk(YvvQ3MM*0I_}$3Q^H@`}Y6Z(W%z@D%)#5{c>88OjSK%=9GwAJ6&3mp_u-<>Nn&->D!t4_2fw@F^uh'
        'UcBAR5BZyu!#U$93?ZSb-2lsk%Ov)Nu?}39jv2RUOKJn+;<HT2pDRRoymmVOB@HJ!Q$9cX*|lAG)l}lLm!nf6^5*lw0?(c>-fLH~>DL9aV<fEDD5gDtw1OnP'
        '$;P4Xg@QwVm`JGK6k!7qQ5@pdCxsZgCMsIY=0#f+?a~2cS`K&V(mTAeEiA@K!d@J7Lva8)myqN;;@%mJ%4@MLx}58}P)c*;n1-{<C>(r(70cnqA%U*Qc;4i`'
        'p;uj%<=o-9VOCvHOPtyvmy)4QZru3WJAr{0JdRdr{E$<*Tv)buImVG4R|51Cah0xIR<@r7;(>lU3%qmIq6*Ksg65fGw{!nlRgS{j`PV=G*YnTPgkgJmjKxK6'
        'K=m~&Hy<_mrZJ{TziBq()FoY)-m_W)Fr)zIsYxXSkT;-Vh&kWAup`VLtPeSbgt~OFxrDG*DRK!4@}~>dYMsa0=NnRV+8Nw8(1sdv!hQpU$$7(+<i{yuwOeHS'
        'BoYe-OjE_#t(_pTgDR(U)(E>2F1JupzKQc>lGF`LRydg@>pV{0S*Zt-(bDRK)3ZQZB#+bq4bO*9QUB~di_>hqe{}XbaGeu>(VXD4LF)6X(`gR9hVIcx<Am>z'
        'P)B{zqh3p(y)KC_ja>^8y6}z#b&HTZMncW8#!gQAy(@<=1y`lUZMd5y?$y#Z?z4l{c<Eg&U{@T}T{6Y_B@hjJV!$kj6w11r#g}RfO7{YZFItQNnU_m8hwbz|'
        '??)5{Lyk-%QIK4^m`UgfB`2nvfrLR;>Ii7!lEhSbS+l%izBvb|mTd2V_n@1qsAHC;Z!6z(C~z-59DLW^rYkR6_{x13#8>yamkJOkh{-6PI}<#g#4@e9E)@&n'
        '64`~1y|$6iSHvzemYoDgvU!vVR}Gm{JrnnT&*o^;7suF43R~zN)-xPp<hl0BRMP*c8#3L#2fHYvIf&E(`Q!KI-KW}{8xM<CZr3%n{BvwqKfF15_48;Sb}Vl0'
        '-hCFc+p3>n)!n(P*V<lUzp7Uu+ag{s2$%t_(9+rsl$qF&=pV{2=IiloG;AP5J2GXhiSG}x@fl*#QyA`!#kI(d_BwUxqb*$7CcB5&>}t&65scx1w!1&ZXxFgb'
        'AHLR~-2C=H`<uoBUhf^eu_V^McN)hSF~GxZU<yFIf`q?NgumN*K+gLwy}jyME=@zwE;A8(oeEf4q*J;Ot1R>!xF*c%>K@qRr&CYfBClUqK~mN(=n5c@ZVBt('
        'mWX&zS6LWPAMA<rAh*CqkEBob`)l+}s_**^A~ahn->JK#<s~GQbKP*mFo^G+vKY(S!>5XGaEeZDiEc=$Zdx{pLrp#)mi3mFj6w`4og2nI)k1dgYP!iZyEUOT'
        '#U`zeL^@U~jUBUOiM;i29-;kkoK$m?*!sn86q@FUlCl!SIxD^uopn~MGKn@+ZEo{rLz^OKNejzT=4rd+TkSqE{T{J;S#y#rH<lLg!hkGmK|l;K2ryPA2$&%X'
        'z4&wq>e6fk8+Dta`#YW1{{yOErw0'
    ),
    'validation': (
        'c-q~4Yjfj9lHhm#iW>MLOnOL&q*^12J~~3~mS=mXE%iy=GjSnn&=64+F#-V^52<ObBJM}rFWg_UnOP51;XzW}j@g@2gxw^NRh5;Mugc0Se)ETeB1;c$!}x&5'
        '4{W{3mq~os`GyU?8?bp2gz?>!75QQ?!6zNZaXLvlU-B$Z{X9v%2R{k}KM#}GU2i(G*1t~YZJhJepXcENXJ!{{0srMo&R$1JoV?k@{wkbjY@LQHKi#l$0CpQC'
        '^ZQQsD(CBroeV~YJ-5@j<m)60p`Q(l{7q8ijDKG9dCr5*p!KiwVOc>U%rhP>SeUU8@Qec<usmTAA{&G$^cQVd0gqGm=@XEJAAI^Wopw45{^%=x!oMxlQ6R&S'
        'H&0Td&2a3kGhPIVm+}Qqc|12?9va=0kMZwhklbcGec%B!Tu0&D&$$=yc~Yz+cm{-h^wYo#(r}TR2;kpnH~@Ozg?RyM>_^^1rf^Wr;v$M_tww`+nq-+V-XbdI'
        'c|q8Dfn2>>vvJ#CTY9VaayZEOYMrFegE?w82l~W!RdC`9SU@0loTo{Y88m%mJOpNdWw*E?NO;E1FFvr8r$x*J$BjS+{5UXP!YgFlY!kyjc{uN=A(`skWp8u*'
        'k`*gxG*4ErFu<D7Z&pNk#ujO^l58|kd@-f*21m>y8DJJ;=R_Aj{k7n^YF6i-DO+G8(a|wsO&=#L=ksL@q>cu-ybR=JUF2D(b2SfnoQI2W&fc8<xgt6O!IaA$'
        '*%33a(M0`^CQ17b{rZ75{5{Cm)a^;-B|RLC$^mN459GhUw2$)8Jl1Gf9cQdR9qQixUhVDg?Z}9kij@!@9akSrfX-rsk|*ROwU+r_qZ0})%Zl3^^(Jr)atlL^'
        '37<htC8{^q)~>gF=|6-?+G$mN**p!`IrFpo4ArmXz>#ZksANWVd3N>j-G?idCQ1H1fsHM&I`cI&IdoY7!Ua4@iX2~%+Amc(OyMVZz#q7+^RwK)<Lq&nWLz?<'
        'q!cQAKLVZ%HlXptEYJF&#A!-$I@eb>EV~caa((VO&<e%1Y+nEk3lRH&Bm`wPj0-LrU)})2Kc->M@oT_wx-5yo4CG5cmv|Yd#DP>6q!I9Kr*o<arE%vAp9CL_'
        '%hg}rzIWMa3`<#prPwxVg6~|O{q6GYhclNsR|w*;6$c>{=-{A$N|K`C61d4V27+ctNYYv{(0%3S^JNDc%!?EjELV_WtZbWQaz1^6wRxW=0VrHr5uO8??lR^B'
        'cK}EJ#$%j~*`ZKAKKyT!g_%%>z?*Viq^wU_kk3YaHXgEH3m-@ivj#T=mhsow&>aqs21EDc<fzY9VLXUuN5j6Q-BB+70mNswz)YyO*R%0(IP3#{JeZBfhaE%K'
        'rwrdISp<9*;T;U*-)2cv0NJx*4Q2&3HINI52SXnnjNIy&fa{v4b8cd*(Q6ucBxdQlM>J6mGmQq1VURD)xsGV6zsF3+#(2O6w}7RX^Z<r~)$(Z#Va!Q&kRN6f'
        'mMvjT1BO!ipr)r;6YT}{qKE{;42_`=jiE_%D71}HPpYHBj1}Bq&cPtXvvD7AIn?70=6<I8%3eyD&c_TTQD{CxM=IY<`n$r1+v2))&Xaf$8aj6E=gaJcy%=U{'
        'FhCLwlpAn>-5un~0oqtkJvG7QD%Kik<y&uJYFrIpBBa_NjW|gDU_fdVM0>DKz#H3;^8EA(byCw?sT~0Tq#p$c`pY^X9s{9&)?lpGIe0*5aHJkXoA&7w=x9_b'
        'xR>aVN%Oo#kryfVJQl7<5d9cbmGA;P9rY}|1IH)h>iax<P#>3mwv57C_4OZF606Uth6fP-x%^Uq4C`;jEucHc38*KVOyCB@0@rTE6=v#tpW*+1OJXkIkidX~'
        '0(%c%#fyBi2CAv&Z^C)r2d5Ai352uH-T_hj>}P+CKz(+x#yR^@pIw2jKx^7(AJCicbUGK8umA1r>ci#fhl@+^uV<H6Z!gZ-j5!bE!7aFB^X0$?znv>q3HTcK'
        'F)=uH*<d^z|7kEh8H^8s=6^kX_x8={2WakHy}o?=9!&L&JrU=eOquh-U4<iL5|ym8_ohNM8eyZo_Y|tt7+ZbmmO^!Mh@I@*Q2_KakpAO`_a8raZ{A*hIjLSw'
        'r!OGV*Hh>P<Po&F`f&Q^Gw<E$k7w^};?RCU4>HhXVRDYFiZoKj#F%`SB2#8tHLH?_#-S+<^k~dLmo*V*%Q#!QvZ$j-!z4MtI|R??GkTM^MI89(3Zs@6Wcr3('
        'f<%XCU!`xCp&Jo%+$1WoFAdosNc~6XF5<$`0tYm_1*MVmXajxm&-182U+<RB{URf6s~j}#>7`;oFa_xQ4(z@&7aYp}_3`cH*&9zRp<wc>De0(dFm_V)TtlNi'
        'b6xl5#$4a7HZ{(0=z-_!3BSrk?d$YwZ?o_Yj1g%0CV+Lgh6%>pTlsOg0K?^icDG<gS~S21@py2@vl<NSMMTu7gxw%`4e-z-FBn8L^H&@tLx@BR_$%GtU6kB{'
        'Fu>Lghp7P`dD0`<g}|S@D9MstP{-cskKT``NLY03k)vub6FnrkCl7=jd9yXp6Hl{%hs5?%&{e|uPg0Yby;2~5*HsQ<;QRI^9(roa&`RL_I7RvOo)*Ldcq-gl'
        '4>l<e4ZPFb>XipV9;O&D@p2e!-*(4J)yhGUXWNzqxK{5eygQz(c%E($VM}`(Ff@aGyV3h|v2H~wRCu*J9aP$<&!e4GX=ymL#`Qr0n;E88s`r1snlTLBHdSN6'
        'k^Fi=-k@ou%BX-hfW9h~$6j3@)g`~)utD=kK5W6_80h=pry=sO7Y4Xi`Ns2s3maf>J>0bE9-E;`4JD*tu?Q*H*vN4w8aT+?VSA+||CdxyVxzrdB|wG}SlP!g'
        '<I+4CD}&o(_zY8MV7W$#*A`C>iXw?_b)5|>!rg-2k*%Mq-D^|Rz=8k3d2O4*@eY;(#uL#EL%+7+IZ-fD5UWtV(xv;C<Rx@5g=F!FlLy$qNh+G+VKhaFSk}Ng'
        'O8Tf*!6L2dYMv}{iFSh?d3qa>Mai1C&iGGye1gjw7=wTI@ikC>5ryj(^D{w$kl@$@XAM5lDgobm3p0*}l3%qJYvb$@2zJibBYC(4c-A1sf?8X(W{>s@uP-jo'
        'z%f62f8@PBJ%586_A}gH2_D7pLf1JOv>Xw6ewOjoZM1QE?7-gaG#w<#ftxxVl;DG}=yWGAOuE<UyuUj8_~wG9eT=*9rl~u}gUVi6pJ~yqSuEvTjlUDTjEFgX'
        '#=c<@zw_rCw(z6q)}P;3(3-TOHPPc3Mi>|&pqaw{&=|>E`(b+P$M*`!mo&D}gr|OXU?MAsr}lOweytj9t2V=zRXb8sMW{GzAmLMzg{bfYDGq$=0T%?xD*R2f'
        ';PC@T5tZG_iA?_46Bz>{#N}K#AXPn8pDJjY80?-jr(((`{3`{qYrY>YPtUL3Qdmy!%wJ!;`}p&D8CP+n5>mbU{uTfA^_$B>^&lVWPb2jy(BPj}Vc<o64jz?$'
        '<e4F?kymzLzZn~Ef=Z*7=VSF;NuGb}fjPSifx6XKw|?yv$nh$PUy1NZC*TVfCyySe=x&b<eqf-2r}A7H);iXKSr?$XDmLpD`Ml@)nfl)CxzGU&3hoL4Wotg4'
        'IT=S;%bcDZ2>i40(VsjFzjQIMJ0&&NXY*wd-vfXB#;14y3(GPZj$g6w=pP%aCA^ejX;*aD?a_<JFke!%*oD=PyH0uw`qa<ZVmZ}@0?$S%+XWBhLU8W8k-xeP'
        '{3$@Wc;?wHG3oW$t>g5j*2&1yE!IeoF7?tg>9lm0{Bsb(rY7%}^F@P%X0AJ51ydb-zD|hSZtSUnDk^ve^TWtc_j8ZX{fwwW;~N5j-Cl{19GK$oE)61j5NPqw'
        'bH2{mX%6yp3)|?kG)>Z4V~07Ssp#d`B7w=7L#r;ugpfuaO24suGzp4rQ9w&hZQ8sjpk#pnl>anyDF3OxO3(uk0Y`D81lDA+6f+7J6n+aca7uIZ*1G`NX8{at'
        'EK}-(f5QHXXw#It3+GH`F0p^Vx;U2zPT3RU_3T)S2=CO=QT)b}$Gp|)>NU&xXJl8JA`u54bgNT@Nh2!R%)#)dKj(7P6*L8PnsGm!FT3fte;4Nl>J;La?7JEK'
        'ckTT0?M?Uk@A%Ijd-UzwK8@UGi^#vrX3)pa);w13U78f@?x=@UTQMF1O`MPCe>uPS+qomBC{B_*bpYiBA<DxokES?3nkP<JuH+JrW<R6_H%a{n6Cj5isQg52'
        'pQp^}xzIjbOTjdj4jS0<v`#_XyUr;Q5KpmWY&CG5o&{(n5l~AQnS3DVv*+^9Yp?<Z=t$5+DbJGVf!o8cK|#?BD-taUw(i&{Eu1$x2O*}ciXsk>XP&h(ljJjw'
        '0}blo0n!`epm32ynB<0HWn^K<qeAEZJ6em(5ils2zk%4akzY|8yO;SU473%Zht2#l6+{t#E2erim0^ovpOSBA<cbKxb|bP39m%l5k`VXfO;<bxQ$df_+H0B7'
        'BUEf!Wr!xdt;VC4VnBKXx(zG{{6l7l{lJE%KEJqXR^&hVAmS3cMHmOrZ3<l=)|eKgEiDX#+{d%2J}pJPebZ;HqK`})`*>6vo9!+3pwy!S0(*(nRrp8f?w^1l'
        '=Slt(ikwj4r2%(dCq*P!3TMfrAjlljCq(1f1<?sPNUNjU{gHtVR+cUoB$4R$TIMG(>X{@><Z;seRIO0bahyvrb4=@@71FF!c2e>@SW^H@amOtZfqdXM++qn%'
        '5m^$HW*KozW?J&5{9eqS5+~VQrZoi#QBWLZ<N2CzCfNzCgj#66L|=<wZ#j>R{Im_vMFv93BmUsW*ykgdWw3!O>I--5D8=L`K?BelQ8*7V-}iCJ(fiT|gZe~$'
        'gAIz<fAGTyN$4t~iF{y9WZ<_PPtdy6)%quLs(7pzsT}Njcw=fs<u=iGbBErgfQ2*whR!IN<r{W7Cz1DKlX@85h(gSx)}1lA6+}nMN2U0{;}HE!SNTomV5QEP'
        'ngp;wh={Y<9?^(KYFT|LVxXOx=QA{Gy>@TS>dvX;^Q`<0Fe1zm@fTTxl9`ZmcH*l0?sYV~0RJMlG=_z2#xW-aOAZZ9^0irnUF>90MXvir%)xdPlr0aMpa2*2'
        'f$PW#gH6Rc6B9>?M@y`@Lh1mY(gl-+><<{okJd|nfV(ovMV2B?+>1Us=J&7sdrGK@MKWkc`Dhyk8zk72Cq*WhX^Db6V-d{OVqq7<z+LCNgWbWigWXZ9o;lX_'
        'npSF!>qL*Xk^4rmd<6Rqgzh09k81W3UxBQ)bEnnC4^|JvvzJAT8{H<Z{TT=nf__p7=Z2(GbUhR9B(QRnJc`6tCDfGRhs)v`){8oQA92jsTLTrTNo7SvL&=?G'
        '{%ol@+#37mXeHHlv(DYUo$_(BpV2mKHW&jmaJv{lS@2(CW*c#OAj~*nZyRYaaPxalxwXU*nPsQdD`O|Ft+dF=K$lHqkI@^sja`LVMio2y1(NNiszh}6O$Mi;'
        'FOV90>I^x$>#xyHc5!J;8-!dcL5SCfebqvd1=b|@ia3NRdts1y^X$P6za&wyiuJz5j>Cn4*q*3d8q8pkWZ5wdG-)-tJIBnvnoFFacQL~pXk#bQ6g-IeBOct$'
        'oO&Q(kcD^kV1lf$K<{vdsf{j1DlX}P4lGz!5J@583~1g7&|QQ)3eYBEJqs*+Wr3QeXNH_okaA)V?5y6DLM%9Ya!eoHJC3Rq@c^F{bYX3WXW&b1>s>0Ex_%Id'
        's8FxE5Xd$NVnvJ2=pR_;=5>=BUm>Jk2J*J@tD7FpC$N||4Z9ziu({s&DNfdKUDm&%vC21Vi~`=pppa^zE}+c%SNFYW^ONQ>SD->}OSWr57B(hE;KQK<1&9Q#'
        'h*XgYMAS#lLSiEpFDpsjj!ZL!rSon&j4jm0MpIcX7iAQQsx3~h^a23{%U^UOBQ9|Lx6{k>8&>M-Fh=8oDYs&UM32bGW-9jghU2g5P7X#R5Z)-#^8g3U3^Irb'
        'o0=85AJ%{e-EGnO0nI>9M``7Koyit!EQ%B9bSuX=xzIbR=EB%Vm9zW4Vw@L_nulEaDSMJno-MO$?FKNdy^H_zN*h%MQz=mpAKmEausY=e^G^P#Y>B6aR%wZ2'
        'YG=D$!!Acq4Q%(R4tX_D&xoQLBtaak5RE;WHF*l{7WuR~Jv@YY)}|EqjeWg|tEW_BEM}FN_@j-j5U0Mq(=n*(&S#qkvo2W{5Why$5Zji$Ab_aYF(&F`CPiLk'
        'U2H|C8hwNS?3C3N5$4@7bQu+!+0iV~S5Jz**p~MVaccC1M58`y&QKK->}mHp6Jg^ng_iol<=@&;Vk0C+h0am5?O3|_s%R&QtJKS>VLybMr7fz25E(HB9J*{4'
        'pzvSdiI^^lNBwNGvwh<xKzm!v#(Z%{;d^x=B3o2A3!}CGNk-kV3E|Yd4+NKyPg_Tqn$sqi@*4cdx_4mFi!iOWomJOrPlHwMS;QVi)CY^~W)%7C>JldmZQ5;M'
        'lGMl&pr>-jAo$OjsuINVjp`KMS%dTp7%2q+B0rh!?HL11<%t1Xr5FHL9~scr3<1>oZom<ZoT;Fdc>b(6FiDc}@<c`ycy)C^P3RG&NlRU#Xb<{ks)}ALO(Gyf'
        'g%SDr{f>Q(B~>JmRCIK$7U8Du#5*8NPDeGAg-aMN!qO(E4o;67z-6sdLksn$V<|YrSqekV3^%<x`zatn?-GK#aX^2)J97{p(MYH|=k>+U@86w$5O4Ib;`r?C'
        'xqd5m(9(7KC3bLgI&;`}Z1SqtF%LG*RCrfcp60@e{^$SszoNqL30L-+zCW2oA3~|BBM&)L!x+X%75q&dWoImyD<Xy@GQtocERfe|%av|WAVnrLS5<6!+6yzG'
        'Y68=1R%dd=3qca4`YS`$JL^?kMzTm|6D9aPxlL>J6eR*b%dW^eP7Qr*Z0cp(>*;8G(|Dl9*ud#Cy(>3jM|xNO1E$nAo_jGY03ikWWVRTL#?NepMJ`oNv#PI='
        'XJq{~&{@<6HWy#pTx3;YV<&bUw(hH;4qJAl3J+0|Y*|q!euh{%mr{;*A8fu}CgGg3>nAzD^No4N=@#oX9y{XlN`Owqe9j1?#toLChAIs%TN-FBR9W&fOblmN'
        'h{mO>4}^t>v04EIn$tbNhf3}|EAsLfsE<LqK_5k!waVU3um1mFsA=r9lAgV^FvN{hj-+Kk^<Uffy}A6~?D8XQ{Uq)(ro-n2sHys{jlQeiT@h+_SNGzC6~AN@'
        'MC}XH?Irg{OqF(4x~+)nJT|0+^x18(Snw1djVPo8jZDm`rWBnYU0~|Ir|zs}Y|NvvOPyf>W)%K>40`sPEDu1B&dd?szCWwM@H91Hu0Fg0pV)@NObtO3>oH6b'
        '4|r#zer1z{e;6yrn1t+BAY&?1OI2Y=4CiNqH)$M%M0y7^r{aDBg<axyf-#nTPU4w&Rp*a>n0HJ45AXjlW5ccCaXhvaytVK{EL9{%00zm$4)Q6TKf53<9Cnb3'
        'J<3E29%op@T&zWu!6V(YYhPtJLgiMoB38n5kPltiu&F+EZ*qsz>MI$ZNJ)v(QrV}^CQCM%uRhwY-rna83&->$q@Q8rcbE|xnOjj_jp5cp89%FFcVh@c5}?t;'
        'ntmP1j~IHRdP^+8+*EJcsN(r#knY>6v+!;2xi8I)isTTd!Q8j@Q}J`48<$;5OKGBf%NKayBF!sC=j(Bv+aintWf3yz#kVtOU&c!HX^83E%QR*j{J7=$#C&dK'
        'C0XEfqXKL&xZDL0te^WC&%H7Mq%Q*xBH&}!y=2BFBrcQoji+TbB;FSxUsM8HbdVy=qI9c(3r;{(tm_q+Q6YYPa--lz!>xd$VF@|f4mm0z$J-&tCFJ3D$itgX'
        '^S1V6%C2pltd|4y!=P|Vb#qd&Q_jGs7{~z7_*R1yx=U;XPHJ!Eat~;^0mxeZI&_Y@Yvj~k#CT~5V&aU4^|zRtEmWDKZlS3I{W^#eH2)6WOVjjz8Qv`gz{yDi'
        '0A1{Y?~jI!2r5RTokiL;*tJ*3SeJ1l*UA<)D;tfB<-u~^eU}wDvo^3*7P4jC<oy$*0wkQN6}EE0Xjpr#q+9hXrCUzCMfQ99^M@+mr9~<)%NeL^Zsf8a?Pz(I'
        'hJn>)<PP^~LN$KGo}&ij%5vGQTx@`CWwi<%MWc1cqobV-YLmo^g!XG?%(G2E@PZ-hOqak98wPkO=Z{{%`J>-Uqm4#it<nY}t6g#~tUVC~|M!ZP_$#y>T}9>5'
        'RaDy*-e{M;`+rs4jbE_n<NBggDc@+YZ8;51H>qp7$#(7-|1RrYQCdQx)s0*25{k`45%NG($7;5jY!S|?T5MEt`--x;REIBEs>3f@DvR`ozYpm}#TUNS<=OZ}'
        'Tq)wG9^TJ~56X*nR2Gzn=N4I3?~}3!X>sJGX!cT`v@NA^g_@AOaWkNGtD0(-$*KN#NUAzfpq&eejjqeHQh(CqF*5>HO=kzyUGo>jeJFp1Qu)kd5RTJ#l+Jd_'
        '1rZxE{3LdU9<i9wlz5gP?E^LJE+I8cHE6lNTFZ7fjpay>+cS#l*+fK6Gn43V^q{){QkF?VGcn&@rI`q?EJf~7CDtr6gm5*i^*39>+wRl-l&Kt{drA+Y_EgX&'
        'Mt$y;X*yJ+fvczPcymwFyJ7{4^`2hJ1%>&v0CY_(SAGsic=O25vRSRc8+oPXzwv1OlL9UmP2Grt7Zy4fI)P-5>#57jOvOd5=54=p9cG{s$yJ0%<w3ft@>@(3'
        'lc79cqT8+H7NY~QS@W3;#)`?zkn7u8+gAs^yuGo5d|~5`rP%_Q184(7K+3)h@|S&<Pcf5t<``Ff(@l$I*xzMxRj0;L0LDE05d9r-RA6QU?V@<iQIlvnAQe%*'
        '+*o&CIBmVOZ(V5S6kt$2<HeP2{n(25Wo>I(x!!y<KHRx{mD?YUMV!+~NH)86{Oq1|9iCqF%R`Bq`sxbM%JTjwuPqiC<GwuF*>UOb&V@`j13V4tGkFOqC~eFj'
        '&QpEXRH?;U0nZ3!1Vj>w2Rf7~YoeHpKHxA}x2I2+G4za&+_o(#`=9mXY$FG(onP3h$gr@5{V*Fjz1n!DQQkd<!Aq5X^=0-VS!-roheG4%SIw1mjb#&nr`@*k'
        'Xn+pisBq7tq!!OJN>ii>H&$7?f|M%Kl%dwE=AqG>-N|Sj8)GBSyq1e-y`42ysESg-JZZR7TbEOI)a+Y(uHmu=@vV7OE{LsXRod_OpG(8V-`nO@8>a!gmoXNl'
        '%GvE(N~78m2kq7s!kK~78r95~v!)8$CKA|<s)OuNuWCbADp$>yFQkLfRF_=n7>6NdZDnDWz%?ZU8x3{w(0kzv$zLF#OD^R#wz^ityu!9mlN5J20l_&P8bP02'
        '*86;;@AZ$t&SmyR{N<(k^_LT~J@{{Zyp11;insLpqAF399`O*%9946HE+|luP%F==YR4czz<O0?i(K^CmMd^fLuZ5;zPSHMbN}M_;Qm--X&$%bXj(B{X5Ts^'
        'lW6ZJqL0ZQeN03j6Qd7tbaQN3qUIVAliCNnu6nUV3g#)yk2Xl+bv${|u#*=LJ82wNc>wY#d#gX8)fLWxw7{`xUDcyM%yh5D)fM;0Rz9RU1wK~UO7x}2{#c32'
        'EGk7eM3>o9eEiIGD4sp3PtU#bu<@p2Or21h)UlnkVLNHWc4A^X*&o}<9@y@eH|Ol!6@sc-MB&{s-<=;viEZSY7=M^;R`F~!xKT-GgI86vPOhsI7G$!8rJs!k'
        '@6AGMGb^uZIqPl8UqGp}At|+5=mdxAyzMMUEviO><H65wudd#n|JkuJ9sshIKBTi84!~O9Mtr3Q(J2pt=UXQ;X_?IA%O*3~doq(Roy??hGIm@RTZt@(W0YeW'
        ')z%|~=C?+<`?t@1qR0eAi5g6HHQ2pA*XI`>&Th03r?MQ_ny9jT%+QD`zpfLOzo_KyrkORsc#i=l<^bOtRW%3FGXSRv`T+mHwPoA9PsL&Z%g-_z{HS+i6`ol*'
        '@$9K_qxUR-<f43AiL_Soao<;enTm#>9~ip)c#c^twHS$0x@s7<%qYw2fyZ>Qxy~zDf%hh+TEr2jn%|JN89rw65p~R!5ngXK#n1Fz_YBqIU)N4<z?dp*KKALs'
        '*hk{cZN9!>O{80;lw8%D8V@NBf+1PTh!?lCRAb;}`hCN#OjTfUDHy_5E*modpRM4lNca-B;>2nftm=$ZG`+0&Fv=9HWt~>H+87J&Hv=0L1x-IcVb&|8+{Zl3'
        'X?aVV@m8GY;AA)G&|uCc{r-hnFClH5q;F~S%+iJWo?e1Cyy;oJmnKmUU|2~5)I99n#sH?xtS)W%8~1KEus3litH^|{T${5OS2PdxY|_`V$BIR5BV7k4Q)~kZ'
        '>)Z@Fse?}1diARHFQD_=x@~}M?^hZ=4|r4BK4@jQtZid#Z`U+lgYQTj;V81j$9-`;goFX7(&)mBnG2$sy59seGXNy0(2?c~vqj5xfNGG$8w1!d8Il!imc)Rz'
        'O^Df`^%#%lKRWb^jA!EZ+cC5gm3Rh4t3Snv>c0U_>a<dUqHb0Pg*&ZC)f~Wn%uoxKS~Cmmdg9bd+=ypMBVFwm6svXC1(n(da^QXF_zP871nk}OjW91gapgn9'
        '!Cdt!3zcMvwLcA6mqv7RCjT`~OG+&;Ys^70<z&{Ja%GCdtT|V~AmYsUZ1kbe=4Wl%++yC@PM&UI=~<nj4aS}A;AWd;XWLm>E|1!<d*Zm6l0R$JuEw1j8lH@D'
        'HQ90-16n8D6d&kVrNwls1ID6;7>|PgfaykuMnF&2jBZQtq*U;ZZlObNW>~E$?&dO8Z)E%1<SJ)Ie^iDeJ@)T)7P;+1i$*5sSY614X@tNLMu*}`bij(HKB6fD'
        'ARfuElF!_7jfL|LCeF_&sA?T6P|+8Ri1+FYb)3}Nd^C2_t(B-cwT(xan0}xyU2&}I+k`o?x-u?Wh)m?h26B_HLvFG^a+4}@CgCdMwif$tiI7xF$X-xHSm}LB'
        'x6bvC?=D{d<?IdC64^*Mg1=VfUvy1fl7e%z#%s2LZ<|gu#2KXd<DRvHWP2D(gQ%Nup6Y8xStV29(%0WAZ=g=024I<w0!u9C@AR9-qg0g=NIOJw&gJEHP0fw`'
        '7U84F(#@ukMAC}@v#U~Cp}p3f)6Nbsj|OJdgmjr()3_-iucB7k-tJylsV<dNCCd=>Z7yiJ3e-@0XY%qZ^M=8@!cPOLGkFLDyh;u{tVIzqm{KNw64$s^&FBk;'
        'F8nAGzuuVPm1@S%%j(1v;%m3F_q;7tS|>?!<<m)jYmtP`>9Yw`szx+{HM_TcNVFMSt0;^?js5-?l~>Wh_*|amZL5ec$FX{(Vfv-VuiJ_gi(aQeGw)EhU#M!U'
        'YPMt4+O@bg6l#GFRF=fo)8ROHCYcbUx+s>GK~a>d7vG!89Dft}rcxQMLm+eu+gB{sf4hDbnK9dxv)&0K$;=6xm3!J-umNOsn+z-CSZh%8imUz=wK|Ta>x2Jw'
        '6}-*;DlLjuQ4qX-nBR_|_~BYNU4h!;-6lEFwe#Zl!m)yWl|1mkym6sS5i)AQE7K4g(<Zi7Ycf+G>b1%Z@i-euPHsj4A9K3Bey7<b-9i&m1@f*?Ah85;tOO&w'
        'iN6r}v0tV|%*tBosim-X-!jwYa#L03E)vx6huw=grj3Ie0QMa#WK;VmMa6Eebom>{_~Vu$UpW_5-+NQPtfxi?lLqzr+OH!(V%@t@`3ZaJ!3SJi6^%u?|0uHZ'
        'MgZ5@fqa!~^a?iaUal|xa>FXCBT}JU$KcaF>`S|C+1sudHU;^-s$uYvumVciRHMww;FvV1!0cxn_BRY)WfXQBgk8qqxwpT!`*H)YXSTp!%e(;b^DdUXuOcPa'
        '&}bx2)XHW%a|(!KPn0B84k2qUdA;dqDBEstIvO^$9qre4)YNvoU)yn0+r$0Z9yYfX7fVl>nuAC^294WQn8|2}aro&+qFK+*DIRt@9avb;!y`7&BkA-!3}<<s'
        '^nitCyxIV_t<FA&$Y_{<-0S=oPTOf?'
    ),
    'positive_controls': (
        'c-q~4+j8T`vEVzu0t2526E;bSx^#OuOP?cmt8ML>>F$s;<8^4XI1mX+aFYNB0Cj1NBK9NpGxjg`PxebTFLeb9AgONocvmsgB!Q}|s;sQMRaVw-zS&t7`Oa06'
        '?!@V>w|ppXvUI=m8?XPnzBkFHNqRl<R^_b!2Yk};eZQ0Clbg6G^Qg@7aG4cJncT+VBumRYn-'
        '_!SLub71uhTh6%Q%lFSk*HDdNcUHyotRx^DNEYKBUninG{~JS}wD^>_aP4XfI3SsaHItP_RrU9odP#=!AD0&68<_<-'
        'AoojTvGdXZdxMCZC~h9#fqxfAHqfL$)dhoz7|XS*Kt9>i|aLosI{8>#=~Pc8BmN7gg|YufJH$%LLlb<04K=>SQ><!IlrpI4s~v946^3&f|0vOE~zKW=RqMdp'
        'L&zX=%PVaKePH{XL%GaKfl4;>Fec!K`%Hhc>U1as^Wy&9QV^5N2zoclSe$Pd$?1us_MOq6negc{0In&*rO1x#H%hV(QE`hx`39Ucf{`gJ#)c0`SvovuF*Yx?'
        'zFWSsurq;}GB_IW#{H5y(|e5KTx_d{4cNLuh_7kHcwJ#6_*(es~?HF_jZE7aBB&V~z7N16<du9zc!7BFcflLN&}ho@BRi{!lA@7%t~o87d++a9k}5^^@5{c$'
        'ej$X7lV0Uw6(X3CwgdOD3LtUME)|*mEF!g(5G4seUMu!h`?be>jI_aJ5S2<-'
        'qf>oV5;N)G(a5@McIKxM(_g9L@XC#mIxjv0Rn#F3%Q{#sT=$rWJa^(^TLJTzC<E0w6)&o;hb78Z>|coRzRfP5Sh(*&;42y!@Viss0(`JS*d+H}rP+%bxslD8'
        'CHVyWzeC@>d)2VhxKy7!4TxUSEIDW0Lcb<hHB2?@4m>PIs5PyU=a;+&kpYL;5|0=W8J$TaBa1jaN+K6lgRX>7|xuQ&@{QVFF$KKk~N&_3co71JtRn8h5J#)-'
        'Q2cAc@a$TA7PJ)B+LXUFXS^sDdBv<7~%k`(CGisJB?oW1y1tC_M1gvh?J#^yIQ!U;WTquC9o`A%bwt7M;$!<R0kFpmgF=AkJKCP^`pgzPyRN%HP=A9qb-'
        '<isX6`g|nzAEARGpJ4`z}f(~FBi#VQ+y}jMNSG(U1|IiWm-2rS~Wzak5$I`?1KrL5EHqWkc8v4cR3dMmp&7-'
        '@bqe~(*)!=Y%zat5xu=`2m*JJAIG4b^T?Gckvj0`3@1EK4c`Kr7@Mw7;XX>W;Q1)gRmy^*x@w6r7rqC;FI3}>rymB(S|B{*p%&<)5)LO0BSe3oA?qr8aK@1m'
        '4nC&jJ&dJ`2lu;}E^zZbCd#P_U_-*N>63di^Q$Er-`>euQDx<82vRjqhXUy-_I$vox(Ajg`_anywzfqK+~j!qLG7g-'
        'Vzw}_2m^t<Qbzn?+yaJgkv;xJ^fkMNUUfaZY9PdtAMfP3Bv>Q`5D0P=lO!0#WTB^K&=AC`#1Xx{VA;$K#X_MUgXg4NXNbUvKE`OERy`RUR5htu$<<I}T~5AR'
        '{_{o8{+=mZf6oJrr{UwvGL@a-y@C*?yQNwUBDZGUecCh_R|kAFFm^nG^n=O2#3cSmRElBz$P9=-'
        'YgI6ON(e(M>;{HG5e&awTE@(oht`%pa=gxxuNbNv43^yI@CJm~sFn|{ypi8{rXLw<r3>Ia>VN9X6qr#P4{t?S(X@7;g>_+P}T_~Xyh-'
        'v=Q+K#bF`KK>G_%^xp@dzXPPD*im^{_ChLZhtNB%U^-Yegz5yM3%oUfBAJ5!=%LduTh?7cfr5-'
        'ouJcs^WhW*`sVcH<9Yc0=!awM^CQ9X3s}2*{>byUR}1}~262EOw$TF3%7<8aBW+ODB1IHyRW(v5@YMO~==+nmN9RB}ErT8ovHs>Nq}k!&KGj{TOHyw>emFZh'
        'Kl!O?0$8vky!9HT>;p4os?(xY%F<GPQ~9m5OEqRGz6iuREmPJ)NktXf6-`;3Wh_>xOIt!@=lH*UJbrV2{8sV>7a8Av+gDQw--'
        '>L08&CU0ynQJm6`R=ownoFZL+3*1Ybss`hN$`Wpb3hg>UB_<Ja6~<GZgfFX?-'
        'fG*$080Mkt;XfZ^*pjo$9}NjvQ`Oa*V5CoqFJmKy|H2cgOGZ~zRnuasg^zx(`J*ZV`AG`SxC?*9UMK;L0PYY5-'
        'E$n|hA1Q=>LI3a)pSO7i&v(s$RPgjepI4^9tKs!h8kG}ug*~wY>_T=<T@nGP3jT}z)XJf7=E3D^avgC%`M-'
        '^{pI8Y3DJvvUl<{!ZOGtLO8RX6b*xI^*PTj+GA@ytuJyD%xTZh+c7&^Vh1c?>k3DurYKrMgll8LY}lFaU*1{O$$=004HzA}Xm+v5Y5Uzlf1s3O^7H6)=u@`1'
        'dSHyOYIqq?C}0j3n`;-2vcLE6|_D-UU4B-UnpI+=IsGLFG#-'
        'QsxiZmJuB$4Q@b<Fob|g#P^eUS$aohnI~7PGCl?=o!2V+p3hYo$}g*|jIE|6NI9Tr@9<TKieDEsv=J$0_AqW>1VRZNaQQCx=q>>r8EXx)Wt?{X{K^l|N}k>5'
        'ZehDu^K9}7+fG2O>&~Oa)ifGGsR7vJ)9!G0@4)+=haW-DyYl^DWSKHIgVhq}uuEV9b82q}H}U;6xdywxD@TMdLsXKh<3Y_5C!*(x$r`x_0g4w)DLD_vJ~<ET'
        'qXGzk6l`TM&tMF*8}3(vv{-'
        '`Okfgv8fH1mHvge_;Y>X<8K=bf^!ik_s?9TlEdIq{G_X{R^N$hguJz~F}d~3{5NQ`tEgF2pr<OQKuB=rGHT^R>023rE}g^mN=56JSMDl`osj_=EEW#A|&{2)'
        '*^c$!Vo;2|)Bt88k}ZhhofxaW1Tt)l`(PLJ9>N>lzBROT6_T9U4`MGYJ*E=FVl&=7DTU8=W4$MEPqS{OzzYCX^$s5D6M5$OE@412Xeu-'
        '2di+bk}N3c!1Y_8uw*m$tnvdP9}zhk6~)i`Z_s-b>M*0p6}*b=Rl>9pAVg&}x)EbZLx_zR)OzI)uiC{Q6{|0-FrlI!{e-x8G#+unfpYQY|wOC1D~RV-vJ(lr'
        '%}6rkH6yP75@vi%Fg=OM~EhYSqZqFtX`iPzo&!{4z4vE3Ja0U^ZA=sPEwE`>cG2jsae;Se{qKb5;^GTJ0=(8lWS(XZUvwPsZ+L9+=aL<-'
        'Qr$F<_y&#&l|vF!PDu3~<(;MhcuK`}kOG;G6vEJ0KSHmpnor;HPf^l`MG)X8y7a<9k@ifRb>UpcJ@@=Od5!F_eiwPcmi<r_`D<e8A%=Rx{8;;oZ#lE*?eIr%'
        'N9ZA%2<Iq~+)N+_AeGm`d+FU&ZDe-vO{BU9L*t!;jRjr$=rDu;2#+=wq_%YK$*wQbbEkKnMh>zOe(?|1OjY27p8;xhg;sT9df6Dj~2Xk|g4)vS7rjSPq1ZVs'
        '<T0QXIzHBxeniAkq|}G<otXQ!<TlMo1NO&xjHxdIFpXGpC}4v<8H&LkuEiKtxk?nW@lPbzp^Q8`kKN8n{@fmgqDw#%NeXpTJ~JZ@cKqMW{}FLCa6vWgL51F@'
        'P^g4$^;;EgwX8;n9Wv@o(pU`SAY7_ka5F-Miz{<F`Iiu|M>6i61_G|EHrj|9<>lygRzS=C_j-m`N$Fcqlgj=H3<1J|3O_1sSwoYlq$-6x1)}5CVx~36St?m-'
        's!70T+7_FX4wMY2pf6kqXrctdA7gTI^POd)_>|mTY)(H&q-'
        'ty&Xx<oknDZ14C70py+~Ge8OqvAH^6~DROm~{xJJQKLs@rjr%{R{{qVH?<*>KV3wKr7Xbd!dxU?b7DehFHp<%GokQCLB@|Xam@V*VSX2GbgNhhezSO}-'
        '0%a{xr^4T})t$TPm?q_J8rW|DZ)4OLL5AE`UKQmuTa{z8@sp3owK8#@n`O?<-+uV<+%APeig9c`hLIP^=XgBSk2QBe8#FHo0|V}f$+Ik#AB;-'
        'U6pgYr6XO}eVmJZQv7~MV<YsV}gS|jS^s9VAyaz--'
        '0=|n9Q2Z6QfR`W!cdLBibd{n-Cspk?VD!`8UIfBm>UAIU$x|=@s*XpBQke3o1L1_rck}E{YNIx5fCZ3eSt-fh<fF|-'
        'qFw<AcL1vdu+V+<v8Vn>TePt52_QV9Pxvi51Wzio|0IYXx&c_iK6yK;crhr8Ti<E&xEDbh-'
        'y!aCW~ys|5R*e2lj3%O$oMNg05hoP&60RNCB<=k!J|<_F#3LpGG?AEu-'
        'e%Fxin^<f%EL&fIejjO7kTmj%G#9mr^h$H(3ITC)N5LQ*xP<^O$EvSkow1M!4%V?}@!z-jy?-'
        'P2&=%q%h`C=#pKrFqhdpnLH3&Wad%vNqAiHn5pM3%2NtMacQHyj%{iCzJJQ>7nm;0UGXVd!iRu0@P!&9_`4WF%>w{+ecp~nx6}%42LH_mkUPic$D}67kqX`D'
        '<h^=pGgxG={@4fO{0|2vUnN{x%-B(A<{yjc5F74@Oi~|Vtx>;Vxp4<2Z#Ef+#C)gDFoKR!4+#6-7zT)_^6R*u$ByCxrH&f%O-'
        'B#_t1VFF%!Rwt3d%0dE_v4Q@BjRd|LM)+=r;Ds8(@?I4Zw}E)#N7Ty<}<i`W8MfPA}!GGZ<Jb*l3o{9~jv3CW7__NCbowlfc*OC`qNPTdVHAp@q_<ybV^sm0'
        '<OGFfXEjsj(6$LX1;CHE>_4aO&4g%Ib2d87B(Ee?2;VfAap%Bj6^LN<pm>Aa|DJ3pHT>{9y_5;P|hnC+Ekwn=1-'
        '9)qJCO63#c+d9Ll(yJ&RtRLtE34d6M5neV($aZLN1vArq=e&7V;Mxv8YD{aG^Vx-'
        '0y&f4j~*M6}@b{mJt_N$H}%N~}tFsayD>lROohXpLbPp0zCd(4A;fnFiQI0g1DpeY?b%i93at5wR6MX%8!FG9KKx~81X;~Hb8(a?7A1@Un0I}-'
        'Z)A}J`=1D!rvIYUU2sdk0fybrWP;L=o#v5)}MlV4uMb^exvhInz;qWiwSBmLjh27i~sV1&+1k7RF-R?|diaw%a4rXUfPLmG+=Y3ms*dLX(aHy-'
        'L@JXi3%=M4q};Sx{C?N@!p^bleDp@^sI;-'
        '_h`)%pk|1$sh?+E6?|D#G~`yXjI{qOOrfXEX`_fQkaZP*guYakt|bAz{OTi2;lXQ1XkEFQBO+%R$sp4+p5L&~gfD{?aDTlhYxJj7*GZFk8*%i>RFNVkVzWPz'
        'AnN1MgG(fT|2~jm3p8X2rh@X!!Vu)Hma|*lU7ER$B-'
        'MhYa{x<8zw%Z&zZw%G3MBUNXgSz!LPRM}YNYutbD8h^EsnJU84NIB5fz2qBg5?Zs%1XAICRYh3JI>Pm99en0kxN`sIOfwK>X2dyspkJUI`xeY5kWd%YoA^?s'
        '3aOEWv&B+l$2j($4^k6sv_Tlg~@<?>k-r~>4Vb9xRf9x$;Q#5Yg!Lg{R*4+(y3vDT71|_rFNaSWj7v)+7LMd^1bOm8Q-'
        'e2rqY8FaWo&;IY3U*W&HdO#TsCK@Q_3|jaj=O+?;V@|5JUN!0N*N}0@)jmzFc|JC#7po7=DC;k=3Z~wTAD2=x|0Dxd#_Y;^r*6NV{)giAcP*xzwBJlGsvN*P'
        '{NT^%E}SU4*)=ia~Ap)MyL#m_52U*C%v(<mV@<hv5$ER&$}Oa_sn@LzGm|@c4fA3tx?5$!A$OEfVM$J9Wza34T0HHQ9_Lnofm=mV1c*`(;^V0#`-'
        'Vq+)Y=DB^!{)H<l3vC>G<c-$P4v<Oe;u-4GR%Br#l2a{JNV6a~-&GNa?J+&2sy8<aXiN5x@u)w7DjjQrq@-bg8Z9u<-}d;{#?2$y=CqpCS6{E=?-'
        'qM?Z4o#Kd?@Rv!w9zCrWMcge4S8?&itYQQ;pjwlccD?%RYJ0<Bh(1YD%C8+`8=~&&U!w10ch@L7z*+a}z^;-'
        'sjHfV1YhVqn$O%<zFez>wP({4S8knzCDAm<NI4zy9_S=WMztwQ}l`<VpA_|0SgSTgwRtjE2KZ|*W!7vo!^{@u?p$-'
        '7Esu`f4%*BS9!KfCmz%EvJ)l+<nRt0jb1<1%JvFJj88iJ_+t-+$pa`F;o^VK3P*iAH#XVS9ZR+H-'
        '+>j71ZWu=PY!9|c`hR_$6zeuT9z*<PH3WW2=G*|S}_I>nIqA%QgbM_MnZ*t||!<wBxc+m{~Z)DP;IrAya?$9zx`-'
        'MJPA$qF4Ml=ASz}lk;l(WO5a0bIZG@O3x10`sev(<CR6N*kf!+UtjUcT8)P<IR9n}PNz(%BDdle9%rBj$id17bXOBvv9Sl&=}V7d4}>*7OYCP>@3^ka<~?P6'
        '(vdwMPOnLfi6=Hg@PY?dZ@T9Z%~Um`v$ZhN&F;{8Zn>g+c&OkH9K%Y3lye*-'
        'yNHq`)WgC7oV3a#_*!i&0%DW8JQy8UU(=_E+ZJoQC&_h8F1d5R_p^h;`0w_6RJK<$MNu5z-'
        '6MCN5MVh2eRwp3{QWYV`~?yL?vj@4brvojq^nGwGuRPYQ`xQarCyPc8Gu`RtqgiInktrgpf6;y%wwYsas-8P-'
        's33&sL^_$n<`a|oCLFl+%rfJYQ17(3JRt^!0=oURrW+1UM@ETI+-95uu;%FFT4R=*+ycvUg|Zcw|ORwf*fgyng5Xe3wS`=!g{WDduBOrE@Hg4<+3E-'
        '4J>w>|H-9C;=hyBH03FXe~?2;S85+1{fwT`sKql4-La!i;1NWh8$zpFc1|?`nTOIob^#&IiQgqqW{4gXxl^)nul+3nEObkcZ5s-'
        'qnLT@#m>CluWxH^XxWG(F@oI7IgbSZd#VlK;x?=v>rzbeY8!fgazr%cqVC@O`w$&ruspQN}WsE_egYl*kjyIxgr{Yp2zo-'
        '`D%);b~+nJrT66NCeq13q%&(M@tg%RoCyzoD&}OGXQYUF{{dz!M!Lo%b^ilMGM+@MBBsviW2_XgqK^5f19*%`mVIo<(7--'
        '{mK*U2Jm~jzf_)%HA}4|c5)l~6(Rv)1n01hc$_u0SU{P5^$5Ro3_Cg^iF=Do=wt#*kF~7FxL6}_?aegbvu*^{v1Q9AC#cJ9XC^#G;O`GSs8KLFWii0LtQAcZ'
        'CcpcihrCVp*{X?iUJq(qK=P{jEs+IdTbZLZ*Iv%x4ufkQz>babpTG73*$X59zHVy@Z<hnCT8exUeg7$oj+0^VS!hi_rGvA@I>fIF?wzSBruY?OMG7M1)U@U>'
        'hHTc^MYLiCh%zqToP{Pw|Xb$AxHucenSg85rYMU0MjwiUPFO;&%`0*HqG_8_6RBM^%pi{@n1{Jo%m17jGjZ~3m2SgTR!Si+>bbPI%uIK_7QLmzvgY9lHE5C<'
        'AQmQp$hZ63HauC?SeyyFoHRHw&R34d_%-?vY>h!ST|6~S+TpT#)3U@#?C#9V$mIjMe0aF4RA)<#lWklN*HicAIsd7UK@(@)F5?ND2(g#-'
        '5TSJmWtf%c%(%8;UJB>s#_Eq?#mgM(gZxu2Q(KG3!;>a;Yxq51^*1);i5DtaXM#6o_K-ZH>+u*SQw9Aae-yA+{vdqZQxb&!+tV!%@k;ARVrOcP^ZA~ZRY93H'
        'COLGCg0*&TV$5a!$NLWUjLC^f9HZ&=gFC>Q>{uAom!HUa_+-M_tLLY5lRbF^@@UbrsZOC)c)naCP0I}~c)xNP=56aduKxoB&md(-pg~_9Uyq;0BrvTP>inv-'
        '+NN7_4B-aw$hoh$DZkNW4-r`cu6iHvT<Db&8VlZ%qdI9w=31xs_Lm1aUHXN~B3)n(%GkE8j$+f`ML<-P0Y~>Q?N_A-Uu=Ji+J95ajRs0ZGPGkX7?RTemlCoI'
        'MOPXhM!&{t(sPA~|*j$q}R5*hj6;{=itX+S1xG^5sZVQUn@qa$~EwUL|2r7)6M+FO+2{ZxjNm5HBxp~B%Jz0h%A3qQlWKdc5+KtW<lA~L>kwT(O(~SM5nLm('
        '9jM3$li7=sqW$q9AhH9&w>7SW7fPP^X=+_De8X(<agVwe@VU6v#9Zxh`Kir`8O-'
        '>~mAcrpe)n|gRxalcI6a4<>;Mbjc3<PJ~Vkd<OWC!DK9P>}5T;6f9pTZNqgAtJ^8Ch-oW#X=(wpHr2@CBNGSif4?wZXzuy<@$UlP5|I7E?Er!mQsw2^dXIc?'
        '29gA?WiLcpKLD42-'
        'w?rrW>WuRp5zBnP096zHMu`UJ6lXvvJWaB#WyV+7WT8#4UX1#h$;9TLXCn~aoynF2^|(PP$;c*>3?)A58h08&6MJURXAI0JJdC$ExlotQhgnk9>6xo(TdD7M'
        '*Vmr<d8Crd2cjXJQu&DN7qV}l)~?E-MtAOBK5>p7{VUq1-!BZYZNOr)5)<z<dgn?6<*ZdZvDZ1bJwpp_&!xtL8i+-|st(qtAFW!>ZT)m?OLj*kDwU3o3rx-a'
        'JM+m+?lE~^*u4~o16hRqd925I6Tu!8q{Mu=x6c62GCMww34PLYa%lA*3tu-'
        'z*XIE093cLQ_ZFoRX&W8PgcnM!R)Vx@{yEzet|nq5^9W3O6Fx*?Oon2BT@i&&=I3TY?7v3>;@MMcXzV_EK4zj2gLZcwk7#Exa_;w6DJ41~+<xeRP$3fi2#$#'
        'f)B&AA%Zj>Z!4uBKA^id}FVEHDj0W~PjjqG>52scdCccp|P#<1HsW3n}TD<?@4S8<6@o#4P|{mm9z}yR`{mv$w1);T9{*a~I!S&$Fwp|2w6q20_4KpL@+cfS'
        '0D9R#auW_`+qZmW`4pUr+3|*nT@On2<f>!kCa0<WS(-'
        ';Itbs8`W$xVm3~lv(s8Gf%S_fi)(c2&Dx!9tep+OZESpPXczCa2YhV|Wy8p}jp40{YFlHLuN&36+M5Qgtp&8N6dKa3(9}vM&m7W$R3X0{)ZwNI*&0Bg;GvsM'
        'aez=|2&?VVwkZx;qhiKJDcD5H(U;Ux>~?PA$*0;3s{yGYQUcLxF@^+U8W|CMB<UoNG3275q-fKj|0`ugynA6S%}|Y-REFUOde@#Fv5}6iMGH|yW=Yvt?7-'
        '|E5PU?D<St~vD#wVIXCbTIvtvgtUHGCU<%E9lFI&XuIV!Z~C@?)D`WKfV0zg&MsdcM^KO^e~BmTp8)T9aJX#JG5*8h`px@v3o`6dx<_NEUUVXxGt>C6vWl(n'
        'eDv#v$yq3{3W+jnACzd%d+`redx8~`O%FFH!`!l~($nn$PAWtvSZ(Fl<=(>U!$%;s5Cst-Jgb~!Q+_G=)*Q|GY2HcoBU8m%iLT1kbj2kYtsB(?C*`w@4)dw='
        '1xS0y}Er^YB5fN?$uQ)OKEY7m?&;L=E*!Dsws+$Glr&0VF@Z9T5@W2$zb)Vwo{oh;{=VE$P-'
        'uA{?v8++K(gVZ^&?dQLMHkGL(bh>yPPZ1x1l{7;Qw!^Wyf;);JNYooZH7mM?Jzk|EBo<-'
        'L&aNrOLngyB{m46Qo~6LrCZAn%CK}td<Tquz;|=X~>F)eDMrULsqYjLO=>b{pIh#ZTVEi9q)5@;`g}J9(1Ee?>1@gYQ6w$1rxA?fNZ|~eG;5pkn#!14P^1-'
        'GRg!v*2>rCY(H}Z~4gYktpmN+PGxH+$Y)fC-IbcK$>H#eQFT&boU)`!o{XcKv=7u)(AD@JBz?lAJHh$sfdQ8Vj=Mr~d6?5KzO$*P{Bb<L7{GRcx7cM-qC<if'
        '==n#5M)BFjoch8UUasIY>Ojf3~NhfOSfY!A@sypan}R%%R%Ro}Z%b!1G)o~y;S&10^1;HtKHzFP)z=mN6AZZp@Wedb#AO}CT<jD07JHT%YL|Bi2*;Q3!nt+O'
        'A+LotZ&<H@SzJobKn1<XiBort8<=E=4`XGg}(%buQDtAN<q7G;V(>!mni2`|L?PU|fye*Gd|Vr0P2mfND=zm9-(`$?KmT;y=@%6|Pxy~3pQR>eG-_VL6ezpT'
        '8S+^o`1>NO_5uN0s=C45>TT}rqgSg+I7oKgm<m`(E~##B;bq3q6nu}IQB9iZ*Qn&}(yo_lt=1g~&f#M1<`dvo4jEamX)AZqC$T^Im!P>2YRA_e~NU>8mMf6u'
        'N^?D)H`7k#1dV}Yw^jxtmFeOO4@DW@b;>1m#W=PyLgJ(A$NFfJrMA&=}_ZIqd<0E({G8NhR_CpaJFY6~o?K6oxwA5^Fsr-A&9M6whcw}l6?-'
        '_KU1p!b8yVh{{z9UXC=m1oZJ=c1h21^L7V=gj_+c{_N<ydA8cw?n~MH_fK^j2S+x%&_TfrFp_hezsSZCJe!B9d5#lf682}@*=tyFAkqM?d+!QqnlQAZDVI{z'
        'mHTPxi-}u5FY88f-fYa!)K7u;d(Lx*$SGb^Gc-'
        'X3X}=tTy!DQ;JH*Z9M;K51spwZ`ZcSt0*zR*7BI^)4IO9Fvr3eGPOBk^hvbFl@8drQ4!TME!n^>%X1{^wDRBqS=VJ_a$+Wz&%ly(OxdlREusf{2x=F5Y_}#b'
        'P)<E~kG2a(5q-xZy4?N&ClDedMs&9x8%+$EL#%8I;RT5~(a%ubao8XJcdN|yktcOEE)-^H{2}WIHzTZgZrr9^pc7iZK+l=yOGuGjLWl8V=`fj2&sD1y7Xc-'
        '4VWFh2me@1s(GA9Ho+^rLkEQ!g@n40z+FvSx2T|FDBMBC**F9q&Dg97)R6u4ibz!&fo{2<~GcG!pVV%|&~7#M@MFd8uEx-G_aCC0jT{er#Q@o1h*9~nRSJ|F'
        '-{QKGrbvo7*DSMuZfe6GcJB#NjM6e-'
        'N1q&vU5F8*Zs=<M?Ve2l1aEFch3tY#8fBm@;J(ye+BJkBkvBKs|5N_15`2ErMvK#ADO+QAdMY3&@=Ne_TZqiZrr&?iHd4W;E|UAZJ-'
        'v;f@*3tfJt1hQ*$$8z_K#IFqrVa<F8h(v3a?^L4L?&JMPC0e-'
        'eseB6eMlOYyS(tCWKY5|Wh{yME&oZZ!=|~;&LczgNeobe_Khn#Q%EI1IplhOYL5_rw(*BKW72a~jwc;=FeECim){A8t^J0K?H1<zdLbj?1=Pm9<b{wg<wzrm'
        ';<SmRux4J^h1F}P>)tHx+NDnurT)VW<GxpyY3EDu+TlcgYp+0JdP$uZCw?{Hg_+&O|S6uq3uACJIV(83aV<URs?^#V*H`<yFdU7H$4XQns=FS+`1+4EyEd!l'
        'mRc$DpLN}m)>%kirUemqS#*A5GQ;&&`q9o<JUNIl9abs|Em*mOXK4%X`t$>o|VtQ*z_bhU>_H;(pOWk8C(umIdp^wR$dmigsV_>|f8I`%t2VAT4Kzf~{f<(+'
        '2rXN<zj&p>_47R*k?a9o{!Z#vrgMbBDCi6L-Gz9a+s||f%Z+a|WNRLU<7vh%-0L^?dNkx88KrD(to0IZ|E{OD^8_0W|<vCC@Kmbb&)GO-'
        '(3#Fl8R^dr$I9%lUVXLk)kG{SiBP2a03^mw5Xp{y+6j})F`ngWV2G(Xr>ls+S=Aw2vhD;NMyXG}KP;|%kC?H^sqj%jL<0DEqd2`!k+k#lFol4GSxcT)VSYA<'
        'pJ`?960{!+`bYK8~E;8KWrLgc&KPFcnOTEl9M&xMdi^t`VWrg%&iJx+;OvTEUA8&SQ)OB7|T=RtC)~G+#(hmn#o_O;@9qY)l!b0=eaJy=Cb?KMx*>J2Z1Kv8'
        '^oA&_uZ;guLusWS8tIBgd1&+~IR-fJf%9TFm%WYO?wKWsYhVE~&<Hr{6$?bS8CI!C<Se4o;pGP$9bCM%*jwqV%AZG>C#<MvlB{S}q7X$&+%uz>PQWo)iHW0c'
        '0MTKxn-;ucf_7d$i3B#{BKyjYnqvEE~`j!$0rr*SLd|3TcE)X(0<(sPiL-x$!fm1$&P*}2+IS>bgs8^~ihp0%RvO%~`@ljZ)MCayMq-uQyo;GC7rgVBV%P}M'
        'yvv4EEp#e~_%@9@3|6zcs(8zXl!o4zwXuC&RNu#z`V`E1ba^uMw!x~(HjaJ56P?d-'
        'vvQBcFdcHrZYE+a#*j+f(P9oRc@<uK0=QM;md|Bz#y2D2oq=m?iERr533_}u$B9>?~!;ouoMO1LeFSBvOd=q^zSZJk}m|~j59TbfUx{r!}y2=@;sL`R;kf7S'
        'E=}3ZX=8JWuMmDR{P1Nc*pkMPB@|N((qwV0Ti%^H^jl8pyKmTwPzB@XrT#qF#qYlZtjMIPg=KEuED_6=HA?94vc0yN5r#i+SiW=VSY9B>J3eGD<U#eKgsvlF'
        '`$tX`bK#k#UQz0{wXoUId8w<)?8=%<X>!wmFMoblb)l`(?Ep9R7Sri@a8@+sTk_M}2IJ=(&d5i3!Jx+(ON~&zXW_khU7D(s~PHQ(kE|X=9Zc7TRUO)$ERyl&'
        ')fUYj2*MWJV@KRs+XB6E7Z1#)-BE*+N`Gt13_?x}sV;*<|Xhi{tE3;YzmerjngYPwnQm0$kVsvHzXve(t=DJa|g3P~7vQ=S6nd&-%A<NA+v=u-'
        '}0&L+WP+l~jbG~A^&hd@Nl<CcYLdt>V;#GNX?Idu%P;6}G*>iHpUX5+^>aeQZ%(w=xmQg}tKqv2|1BC*Y`e9R|3F}@uXJ#ynUM`{;+H#n1Gn>_waAPZh8(xY'
        '1V43n2Z~LWwJPkbG*^2)oz>d&eZxIgVp2PSYy0#4Fj)ORNIBy*1CN*4oHP8kDNKA#MNoJwXDG++z<CAfGfG?@EU>LaG0D>`R1{9pKDrkGkM6P7H+Y{sMgH<q'
        'r_bVmVu7L`IrcAswjNv^hVvkAYM|LH4t2eGPOF0Ll-'
        'd0jUTR&CY#|~zd)mdnhlrqVBx7i%XVjk<9wY&+&H@75sl?CW`4v_slop6?u%@%kqmOupOhvkbR<?^R{Nxfo)NtufZ#7Z(+k*gNKWk8Pk9WE=KTv*<l=Xvjtq'
        'j=F|5-'
        '&?6+zM2XwYY}_N@|}eiY+=kJqVRQD9|S9YXUO1{wr`XL7yigN%IH870oxlajIbi7Ga3#Tm~T+d52j9R=}Exhp~2zDX8p`1N(AtF^3_I4s>M4bqYP#r2}@FaV'
        'rHgjiRLdy<@sKQKo-'
        'p^PZN3cmPraloWhobtG^C57*}Yt>)`=_o~i(^(~P`^h`}=8ZUjvy>jY&>AdHK;>c3SJju<<o%vN<U2*3~OlHo0YhMa^bGquF8{PlQMu~qJIC!)3mODKvcN$3'
        'A_eFZf#Wy&YgsuEC?x9QHMsh#)DXd(^V24$t;FKlvq<je7l8W|`%f)ik*p5hQTA5oi`Ji5~c4o~Pj!LJIUHw~QbkFiQ{v7+2+ka~^SOoj@s|m+?TV37K0(!6'
        'i%cKsidZ9}z0FjTZx5+gxO@P?%z=2D)JF20Ej`|ZZn0EWE*TuA9T833M+|IHYko^WA`x^quR;7U8Al$XJrE%3)wjTz3?fRQnH0?pGJ3QE?a&*U`_7qe0dW)L'
        'UZZxXD+M>F&C3ST$CR^%Xwy-61PP8|Cfw$uxouB;a<olEJzg>DiN(Bh7e?n!3#oKpMJ3wErDfRt!@J>B~)IcfGgr<<-'
        'omI`if@Yp*_}*IkqM2U%To@iZt#mi?a`xV}v(PHDl%BoTT5)>Ydj>ih6osBs^=Opm@^70vDq08QknT>d&<ZJN!ADhT+E`&~StPZ^ajd9W*$H<bvsn`dWAZy?'
        '0O3V!4fd#Z3ShTD*-tUTkZz-fwQgN9vWaT5P5;H|&zeTva1}^WG>u`EiXINkJ3vVBHp#2Y_7(CBVJtPkv+N!=kpgd;iF1H*UrA)#6}jNst#ga|YDva<KTRa%'
        'Z$RF628%L~VrlL)Yv0LPb1kt(h1a<r#cfC%O#7xzmaTj&np+zb(J*BqbH8@Pj2cQh=o-'
        'B*7fzSPWzYbFXRE1!f%j5go{IaCLoPws+?QSBnnf;wTSaJ|kN3=b8E8WlJNa|!2qYDpYwfM{z8zeMx2D&F3Duk$ncUakWUpfeg^y4Se&jvkgD2lEUWvd^xPs'
        '!QpR!iX72v?t>h8>H<-e*pVEYfdsHeMHOS5Zr15qA>zA*CMf2b7IN72>hgls@(?AF^9)s|>S#zmV@pLS`g74$Yuy@KAgsn+6cn;m_92JXs;TE-'
        '?$wVK+MarDmErFruHjR8K|yF5<~G}N);G4coWsf?=kwO_nF{hN3C<NHex6f$Cf?o8DlG0tw}@s7ZEc%PlYeU_9WmMu^aH-'
        '7C#wO#uX59rE6Y^6ZJXl!CuxV+9VUT{a=742q*vw3#+uinMqcG64WX~6(PftK9>G`%}3I2&{F0d#R^-lLCf?t%VcdxQT;=-'
        '>bSzx*d>;`!aq1_$i*6}B?cJ>;1BN?{pHugh=e1!DTD7v$%@m7t)I5rfsp2DQ1zF<C}Vsd{1atpLyYmK<=!y_+JupEDiRN{JW78TRNp_J+hSS^YMIoH}GHuI'
        'iGSg^uXv2-0HEm~zPJ_%i)n<a&CZON*mmjq7hhj?LAAzD&w_JoX)f)I)GuV1~uQ9U8m(*i%CutFd=%^{E(N%={l19Tu+?C<0D-'
        'i71d5BiC><I*d(Aj>;X=idG2Vhi0HWo$es7k=S?;Sfdrm;m0a+`AKmsf!VQs=g#I)^IYo6&IXLwob3Zx+`hE0=sVlmSN3e**m_^sSMYs3%h&ZR-'
        '<IXeYV=`k?6ca+N42?6YD*u~Ha;iqV`@5TXz*d&iY8Ka`%lc1hi#`1DNBC~oNC-'
        'u)w`mD%$~)#>*7`>gL@Q5%0}K$o+=xy9kM)iHp1O}ZqBz7bhGIP2D5zIQwIYrr-'
        '?eOHZm+v#cc!QJg9vt&zc}>TviadY`z+yYFU~yE5{Y9Xu9Vh2Sk6K{J<8(dnBWTnqkgm<2+UKf%9(?tNe8viAwXfhY=@wL$4x0P(W1x%kW#0cUrZoNVQ6(QB'
        '_dto5kA2rJ4=1vyIUFv-'
        '6|xk4J7PFPu{<Q8&sgIva#y>sdO6HGh~dvR$dD;Ka};y0ygg=5T|$njt>7)(>oHPiL8U>(+@~2dx9Cr+j1L4DX&T<oP3|E%G2o80@Uwg61qpSD9}S;R|UjVt'
        'u^0vLZ}^x}Y>hzJb8hVwC~bg*6(=eTrw~WANHK4(b{OC1dQLa*Cqf?qKf_{yE(34F|gpA#5D3H>~G6SUhAK$<!VG-d-'
        'J2na!byxSMJh8rwbm0_iDJU}8SjX9-'
        'S~Z@&tnDl_FQZYDDuc=0BZa?g;I&BOZ=At~tRpF}=>U>3QBXxt~(DtWc@^su2T_9zkUG@c{{iQqU_hWWr@rn?V$Xp=t_xlmpCTm$p+YtBYZoXfmu5hZEmqD9'
        'U^<I(!<d6X%%$L^&hHCr|Fv64Tb;$YZeNVPn>CSkfR$Dl)j<fGx!i(5sKSQ<iAML@rs#ltxm$Z9933NWoVICKZAYUQrg=V%XFgQ)I!HK-'
        '59GTaf;eZB6mnf{=4h?(!EbeXqSBpa!4SD{K+C|ajTSQ<pD7O{~gq2(+yj=LzwbJpYfj?2tHr7}3NVmNqpJ#?=ir9q?9TFV?6?NMYjMH*s)YOP_^_a%JJ&A4'
        ';0bGXAnjmpEFfYCg~t+J@pEWw6gAjdvf*Ptj>iKxLMsHsW3iu9PPGi`cv>-^<{3SDlZAZ?-~xe~RyM!!Z?YEqXrP?;nmT#+-'
        '|DNgnHn%e0lsp1HfQxy5CD6@s{E@!G)nAlEabEfCjJTxkq*A)RaiW8e}C~jz5xHEu1n@YC3j_iP!jZ~_cLVl>YjIB(y=Rw-t532yGs+6AHx-zx3+|qgI+tsR'
        'ltyVdQ_$3Mz(bPuj)N@s;R<){;(VR&xqOoUaSfqOh1hi;d2PA`xwA1onj3;j#SgT>_n~gTthmgCt61zWJ^J%$f$tENgoybhnhpx<0eIIDWNaYK%l3Z`D64e<'
        'vMdk}NkQA<f6(w$RI}M5J{AZ{}B5K_<==Eq#p=KOO`AUjXjYuVJYHBs8-'
        'K{tA$A5yt^bgQ199iQ{gqZ3@b+Nvja^cJ@|4c>0ojKkLR>Wb9D;Ym^;GvQEd+gQ`47wh->Kj7Pjm4(}lE|YnR*6?Tj+n-'
        '=Gp<p^;EYX~;X0W>(tWJNe8&~NImUnNS3%{F$;m=w$5yE^=gBn&2+|G+XTlh`l)lXmq0cZRaH@7x?8wIrGSNRv&EfV#Zktb|K2zg42wPxne!%~W1=ku4=)lx'
        '~7}o#y3Q*g%8QZHuha5~sXIMp6Eg1!4MQ=3k(YWtKGoO+p#saqOuFYY%olgdzpsM3gOE`UO7qk&krYanC307E~9p7wM$XWF(?F}7nxHp7)ER1x4L3{dqV>R-'
        'LT$M$O7PIA65t9`&26gn0otir+A$Iu2K*Ev^P(uwT<y<q)9VI6|*Y7JG{$ux*4xhQN#3XKoOvwl&-&zYVw0UGL*=KDbZIOjWroQqmA-'
        '3&ai4al>S}1wf=DQ__VIz>u<=iv*;wF11hsNGX+i<@`)|s99XW2vf$y^6p?4-'
        'P6fkyFy$!BfRHw&=ZeU{<@i7j~p<hBy)l}lnBu9sNrq*1j$n+c+7t*<15UTrFaIC5DT%W{xA7PI?m!lz0C^1dTQ3w6m-'
        'MbLVo^6DQeR9=0NP*JV}IZodq)fW;davi;RID(pLSLhzmUVE4UZbr*|V74316|x|1Cg8!<3Zz+o$ya3X_vuW1FK(~DATB@+QeUyE>3LIW_bO}wQ_JwL(Yl$)'
        'YviPx2)b8>plcuNmxww!ZT}=v?o0%{Zz14bvtLkO(CaPdpg=RLpuH2$1>15_Zskg1@qpx=_`bEgd+ms;du7GdRUvOCxvI6klH_{5q2#JE&aWq=8m4i_n8Lc4'
        's$yfkn0oz>6;rRjNKBd75JNTbG4O+}Wss!p7neaed#-'
        'SErkj4d&`?rEbhjRa#;J<>XhYBWdTdLL*;bffZ$q8v%{Vc*u_$>RHUg;%62{GDGNh3^Z6ZTn8#1JAh+iT><P`mrND#h5eG>t){*>v(+Yj%LFRM3MY<<wQeqy'
        'uzq^4!+i9@00qaE+JmN6}rSUiIJ+bXcr3daUul6M1~-'
        'Cgsxc?oDz+}7Z{6Nk)<(~@H5D<>lLIa6yFnbRxQ4x8F|*T(=mFzaE6v>S%R*1`~T`iG(5YdrU74?s$b<NKuSa>oCl^M3*S<tE?'
    ),
    'controls': (
        'c-rke?Q-L|k^l7+2y;~{?n<=fuVl1QswQJ^;@cT}&h}=1cpVNULKZU=$q>|dygpamL);VIlU#S>g8)Bdd$U`~rEaQbMkE@IMx(pY-DvRFzuMV|a%Y}oJ3RZu'
        'iu-DnXZwR+v+<W>w#efoyPdL4wH*Hz-wZs@8|39;#YI&{RbGaRJgdq)6+v-7nDzV_3{JC(m(ik1K5=GuXG{27tvFjIS(LIjX`bcB_gS<~7J@0%JA!8-fv!n)'
        '&*o{q_&69|R=f~waxk78`N3dt^10xPipQ+vMGmledC$`5KHpS>aqpkO)v7^EB3J^1W-(x@a+V^saZ&>6^qy_tb;&+_0G{xj4<Dw}!GOUZdZb_Yr<XMvV0Qwq'
        'ORLB3URVgeiSw}JOJ4G9VK?rFpQ1F0qbkWW>-B-bUBARXz!02337tY>1au35oXw)zl-#C+pby6}uV8l3JV}AM;XH>)M&*5IalF;(XdJ;(-T%nLd@gV(8l+m4'
        'T&(gm4$FML5mm+o&6}94h?R?nPhzv?f<<Y`qxgQj$k#>6art2RZ<A^hqcH3cpX%ukV1P|piGZC}EaRVGF-l%;GQpHEj3@gn%3{^F5ikNumgFu2B$42iSx=Cf'
        '%SFPoDp@8A7V`p^JzLzf+cGLvTf}B8*Bu9ClL`Au{n~f9^q<btr}{|(AG`$kkg?yvu=44UfpBGVWrZI6m7(G%(sQJS5YZY`vCMzu8QT;fI6PVd#>J)t!K?^)'
        '4D=+pUTb(G7iOR)h58e*MU--sqOlf0w!^eg(o*Z!S}JhOl5EAxq@ooY$YEeH=O~V#5DZ0u;IycW7qFO2Z03K3#RSy^O)UivG%1HjiZm;qc|4BNVimC={WseU'
        'b`O1)7hsQ&46r(u9FOCC9VOYoHaR5jA`u4r5&IOC35ZvofnC61Y-}+3B8LexqY^9!>IS>PWOu}HbozJikR}AbCM3~u4hXZn0-;_e6<DtdIc``GxHfdoSJ9^='
        'F9*GbO)Sc!00iQrKqITH6^zxRZMCAYi<8Um-d<g@GS928soP!v8X_tHBmit6%dHk3tS|#xTWnwfhcRD9pr=z}3<(s1qm@2gf_R83y5$O5DI;2IWeXVpoQO({'
        's57WA*<_dlkJSXOFK?Lmm=sEUKXM=dVQWmAtjafFg<~Z!pyAnu%gz@!z}`C$G>(rkM?UCPjH+JRvRIlLmQ;hm8^u4Z6=M`V<cj8iVgT+Sm*1Yght<tA&tZSS'
        'K~g=47+Am;C*NP3UY!7=-j{g{CL580$pZRsMb|=Wy_?eHlnTMTwf>L+6`P{SOSO*cB-2VG7|1L$$~a^DHHcPM(=6O}5XcqGB(4}qI<_9vt?Q^#>?bHQo<Xs<'
        'ySs~{uQ$?w2<x?0={d*@?xTfj6dW+l*ko_NR&qISsoXQtEC~yec^8#2(g!_i!f=*{yB(aEvAv&0!jVOSAQ~DvPX8$q%>Mq5AXyEAs-FC4WaLl_J7^gSN}Axh'
        'Nz-YO3ns$S2bleL0u)7xq-~b8Po$hk=^ORBM1wQ1x5-cJ&3-U@*m$fwgE85qog%7MJ5|1;cKpybbD;uty4$AhS-p_~99+Tin9LRECv0|U0?q@PrVk%b&AVEa'
        '?EwH#cjPzV69dqMv9vrz33dWKoMUAbU|BwV0Q2|ZLqMMQ0QTlJ3&Z6G{kbs2ofwoGR1J0xsY`=_eplWW5p0V3c_BXOpCGKOG@0wizY5sI)$h{40|@_IKWwTb'
        'HBXy4kh(zP^_zPkX#=+LDgmRf8uT+*bolqjJmV4$*$gCTu=nsN8>;&PIHunpCyQzX+at)xJOw0gfuSSzE-Dac#Lf%3>5kYX|Jw$Ib;PbVU}^`0!TH6TKb>4&'
        'UA(?JzX<<)a&dWjeg?YWeVU9(4aQoU#(TSaN8{aB<Gp?0@SF3qtBdou;pLl))Az6m&e((8IS;1Hd(yx5MiTfCfuHGK>pn-==SAK%^n8Sges}f$yQ}c{^x{R-'
        'zL?YoolnsxiMzad{l}B=?dw0Bylqftd_uDa;XbP@S(l?|6Kwl++g&-<GlNB$3o*8Qv~Db{IiZV_|M>3o;^a7#Gn34}oRXwo!<>`=UIU*;EC_;|8)w=>V|TpW'
        '-4J{VG-}J_Hb5m*^x3T<xdqFF!qowwJ<|~4sWq%t06N|PDUyPxNyfu<lqE|rbWNDaK8OS_CWq3xf+=o69E3)Uls6I5mS6$;_1gT?%@BfUwY`Q?B*Hlec#^f?'
        'eiJTJ)D0=q4r~~=*)j<$j#diT0@Aw9Gm7jx_-nuvN|A!F_00&lZQ!5c_R&8hw-1`!CQS=TSNk~s+b+(x`*p8^4foO&p`{QuiWc;5%rYphATdo)Um3fT_lMz|'
        '*JsD4$FHwW(BR4S!R>A69gZDmq7LN4NC>{3r}v%@!Lru_1r|RlVVA4am|aN6HTd|&NPX~Y-ouW0^7@^-o4f*h`AQjl4mO^Hqqu;LJ&gGxM{5_?YhETxP2gmL'
        'JH=!yxm)*g)><X|beDbXvOhKpWW$ODIIr2Qw8Mr{ZIGX3f?xft8Z$A|M&ab~Ef+1(-V2K|A+a&A-8yg>s=cot(69OW^UlIza{o*YY65NeKFzcyWVVlX#|~4V'
        'dm&n}SuPVw(=~+OujS=uWW`~_GV~hhiQF^~x;9N|3ox>0S$-ELA|Lu}{Cf!e%c%+kVf4(n#sNSLwRaqBs)Zj!LO&1v01$vgKvYdYqTq{}CqVJCSa`k~NUWm0'
        '!=rGSq<o03^^~l`h%Hu|>?6$hM?S@<9%I4Xy#w|o{pV7(f<_e>1QOlQr-r+vT9LmV07tT+SI)trK#;gxO$|pd86>m2m_C5j9;VTH9!FDv3NW)XtSRwFZ0>pf'
        ')JaLMf=z)884?!X=F=)z@y~GrB3TXV`SLUk|C%S+aIuc3#{Iv}Np^1Bse$NW#Eb@E_}M37KO>>g^@aed^1dcV4NLJ)Kto`ENfQ6O;6=q=S5=wJHx)lA%e-uL'
        '^_V9+N*sUN<S;RNV2Myp1-TQVBw7nc<dZog8>S{T?OqOV#9-=lQZw`Dq^7+-fB_hQqXcOVtFc_k`6NrqIVK{>AOr;l4gqomvHqg67O0F8h<E=Sr5icn;nFiF'
        'F8_WBVU;9s${vueN6%S4Xw!>EDalaf%bL?Qt5xyOD5x|=Vh+9-Hm3-aMpDGg!;K|ca5d^0x`Jv89+ivLu>9gLGLX{|;1A-<8T`HWe)!^Mc>Ndr?<=34z8KNS'
        'BeqPVTQLJ1@0@w8gWD1|yWzw~uB{nQfhXQ~XWyQke}CqwDe8zn%s_V`mcMVRBu#am=jpUp@Mz-Tgi%kQq9S?)3j_)P|B)84-?Lp?B+oCqW$_yjHEzt2A}UVG'
        'ewrF&AgrV!%I=5qEhsp>)-o&;7onLwHm1>>r&C;3!jA(f0baqx*j2gVwuIfmx+Ph$sX(|MsP7}V0ytd|8fA$riG>cX0fAjMtA^fdq8k$@C76Vy<Z<Bn4xp0}'
        'rk*e|`BH$>5BTQ-=8CCj;57t?dHRXF!xx4tE(Mxz1DB+wcWg4r8Rg9(v_@5V^cw85@=$da-7HYeK$ls^R`zNMkt~;os`_+AC%mp76w1e&7GJw~sy8Xwslkb7'
        'vy9j!bNWh%aG)*@N*i;QSFa7bq^&6&Ei4qv7AwB^80s}9T@gdf|LL7iH(7fFOj&LW(mh0pml-8$QRa1>2tjAwrU7@38=a*EyeV36gT5)4Jr+QwCb4%zYSef^'
        'RNx8Z5CQD@reNX~h|gv-wmZ^2%=SV2wq5$6TJx~xM?-?Fi8N!DsiPM*CP}_M-#%;CqZmduQ1%wp3{8NX>?)y@_O&%Ol2}Mi*CYhDX+9r%U+z2|O6@!u@c7S~'
        '3IW|{VaDY^VreZVEp5T3>)o6Fi8?5L?_ePLzu08x<VneQnCQI&E+ingGIejC1G5T~qe5Di0G3(E3_iMwbh=$1`DtrF%=F9vO5})yZ(!b_hTkBv1|8m;e3MEQ'
        '1}_b~HBHM?yk=m|^E*}#2sYkc*GJjnH)L8LL6fa$3-ySh6Hbn8L%usi*NJSevE$94Q#hMM;4#Dy|7ogI9d6<)Wj=IL8A%2rOg{|4t+EZzw!3zpmsm#23)3j-'
        'QySP|Y`oejlRlW3HmS!5`O(-0WY36Ycz6pM*hq=v_P-W3g$2#jn)nX&ZJNfS_G%giqDSzGT16Gza^cWBY-}{zX_sJX3Duc)El7QqvB%Zc-KK2-$*;B#XIE#;'
        'o;1$&{Pa-Qo(1dSFxK;S=uoRUbwitD`-dm=AEns^zMDQ{=rcU7dBK>{*WBKLymo@dXbcGP2aedGzztuG*kl5K_bow*Eq?YI25Ay%wrcXsu8ES`<vPS8*2*YJ'
        'Ve43<rH3sK_|GGfK+}nzfO9ydAD-x=h^8)P+qcAB_tF6to+5@ACv3#f*v1ha@&c8{Xp)XRE&bAMgOQ!ELE<vInFN}?oZ%ABh;^T**sRAuYx@8scO#YpDZ}mf'
        'q>Yw^8#)k0QFJ3y`KS%cA}?bgTQyG2p`&&Vjg}~gkP09yICL}}z{T=5q%dWsmWFnRcJ6Ah-8s={eB1ooyQURjNN;qG!Sr#UU<_yoZaG9yE-ySEV=1)WNWF`@'
        'i#<xOd$IW~3=Y?b)+WYTbHxmfrfqZSYhqTRFl{GPgfpY4cWA0}Fg|fR!Q82rv&lTv4v`dpZ>sB->K3E!(^_QfR;}1mEyZ?Ece04UePh8knFYYjR~EE0Hvn~>'
        'EjZ$iH<etH?;lOvmFI%jk$yyFN<+m|ZE@()u{xUTV40@+{NAxS5<(VAP}X2HaMMD{iJ%peh^{sJZ~S%&4s^*Im>5*nUw|ZeSfSyy!Rd7exGLG`=%HJ58L*O}'
        'nM+}9_5HTsX&}k(8iLSxV8Gs<cu0@zL<y&%_vZZF`?n`ovc(K5vr(sKrd1--<z<M@?Vw|E3%7TnCZlpGKV%;Jig`1SeaU`%;1BGa#)%RDQ!Pj31K^4N&wu^5'
        'EDMxd7(GAORRGcyYkyP?)iH*ZHlaULPj6%n)FtuI=sZFZ2OXd{;D&2yw}3W68mXqn`A+9^IYZsO*}0{HK9>J@o9ttj-({Zfdis+4wKUb680}s*Vr?Tv?E1!c'
        'f>oWP=t_XSQwVi6;<S<evO<4_=6X8WyJ=f?g^#hZH)19&VNv0vCLosG>%*}GBta06$UWG-#*@8Aw$4k;-MWJ3lzE~Er9Aq?A*#hfHL8{BF7tGy%BlB6EN;c8'
        '!KGx_)OqRzwv%Ujv&xeN@5`eFo1(z{7p}bmu#%mKiSq2Qd4M90P}Ls0v!6H)SXg$6m@c}87_{sbt_I<dYA%(zwofvgCEsGbhDmlGnSlS4-<oryyIB8>%X@T8'
        'rew0nB-fsRS+9R%Xe!)O8CQ$00sD<dgBX?2s2jGG9sv%s2YB+g{wVX+_#h~mGLrjkn18f#D|9e|-Q>B*UB|h|h|M?4B`@*Sq;*JwHErm>nCiCRJ)v^My7Lfe'
        '<oaB9@r<0HyV!SIL=}Uq&Fm3Q-=DN#AV#ubF0YQkjdh{$s6w*Bc@5LU!-3gk)X42nnW_-HtZ15nq{Q~^oU=SE@xaMhuwB*$!H~=G%xjvNRY2lLAazr@AVIH&'
        'g!nE>s-bVh0UCcbW4l{(VxWgoh}kj%&5GIZp<F!r0f;l7f_2Gn>$LiOlccd0dNFkKAKC}0Qy`;3BY<~Z7hb$-zQBX{#+!Znjg`>XuN#>oN-`@2R5n;1EMq`v'
        'm!&eR%I|t|a8rXGZU;TIK##VA9^DMO-CQhOyLN^2UIhFzUu-H_m+6f+-pC|rfIe3&z42J`YmCwINHMP)X5WU{Z^1|xK(`;X+Sf@o)>YDDEUndLglKj5y-xU6'
        'DWcc{EX&%*aW$63?+AL(4XP^Pbr(%m5AqDjU1KBGba%Cz^6q{htoX*($(}uewH#QY_oB5tY_FxACmPur)48IrtOb4bi^+u)hw)a~aAoJP)lOAi$UQ8~Z}V2w'
        '<-pPthr8FD<m_<Q+bWo}e!tiXAH8sekN$O5ST5{XPVEUPR<P!t!+5J~xUzH9EjvOA&t$h2BKm(%|L72!di?+U;-|9&wDp@3QKurQ%<qKde0HS*ZW><ZGT~|L'
        'o)1!z6v?1U0sJ8&C;TN-68KB{^^~Etvm07Cn#n0dmt(x;t-<Yrz(DnqJ1A8?`~;8GD%C-kG~;*Z5`n|?jItf#Bb_6OPXT%j-_wiX5*SHb!DALNlMM13-J}Dy'
        ';!$j(1``10cXYIhYvC-g-4j}e4<DxpMFHtYV2mdnUn=vGqdvOc7F{_hWx<=iffBqfZ#Q^d^gT7!M|3>^Lz`CELxGfVk%IR+Yjsc;WyGI&T708{>$L`Fhn5&t'
        'cjld&t3<r}qPe+Aj29XMUCk${4Wc*>BQ@&KbMG9XiM%iPOa)AmrWp!wUvJ9`dtv8+rHhUGMV62+C|DJ<YcCQXNz^e;9Mw_bjm;!Eqy*hMX$aWH7$1+ARI(YU'
        '9fa>;RhrPdQKx~G+)^6AS5mj-%3D-sS6%}8frE*mU*D{kXQ!&zB&4)=7+OUpN+E66-a007k8XXF4#ZIc=X@CAAUD75w9YyUoIE57Z*o0-qgG+_x(+lNBeiJj'
        'F0W}-7T9$>cJH2ArM@#|_jJ->m>oH^_+>ZFx(b@w$DS@5Tm;PlpHY_OLN_Xo?8X;z#OkUTVk3LzP5}Fgst9o7FWhj1hd<oP00sXX^v+2&46yys+ka+;*Xyy6'
        '1-L6|a=WUY&fi{GwI-A`ZHs`;^jbV2NQdXmJTBo=SzJrGot!XzDWT2?d*|O$lI+PdY?(THcY1kwdiF=Jjr7%4q~!m$(Y3PjHpa>1DU5i<<z*=MpG3L;1WLmB'
        'Ykj|8SG!qXJaBIk<1A3s^zbVA9JCH@w#>|bK`oS~Ye5%ZF}KO)d=YH~H{{T_Tr~Et>#{!2!34rT(#6^>+RW*@-$<rynbH0)HKYB%e@6S=Ga9t(R&+;s5Mc*{'
        'ew_`|+cKs)n>s*65x#F_%KD*!(=G7*E#UhGTuJF4(34;~&nIQ{ZpAaZM2?Qz<+kv(5$mK9JYCkG-pKH}d{P;)wVBd+3U1XMs`Gq9wAI(|+V9$38wG1#7%Dpy'
        'UR}ICyF9%*JwLk)!}`EcPIuOs=d>K6nfia!Rird)&`|@?Le3f?y<D$ZGm4k%GHV9zmKmQ}zqOC43+9%|t+}8@%3xpm_ohwf4`ZHDH?VaK-gW)kK4b@h0*@_>'
        'Du6ZDbexm*dNt@Myx|*JptdtFA?^LksR>ZKl>)1GFAGRZJHR}+Hp13JiNF@JEdnboJpO8g^R=={eDv&h04W`<p(c|~u>B_3z6JK!3c^KelDebcY@l8!p$a<{'
        'v5eXvwFGorH|gpi?>ZrT;|^5=__1}2Dyt;reFIZ^sH-D;X05+CWiF1jms#B`D>=kZNybHJHkJ~LGGwvS6J7K+5#o1}LX-E|()+Of{+PVo3c&OslMdLUP-VNO'
        '?7HWZCpX>w=;KBrlwUUr&DV{#5w;h?CP?P3P~VG*()7OF<58!e@ZVv>&EDNLrMF@8lOm?v%9nW$dEdz3rAt%sW?hJ(1RlX^$3XsUZ^Q)N3;xLOrEl)z?)!s%'
        'rLc7XX`4Qo3Zo;SRtcVU>%fGUp{i4`i+6QxnJz|k$|xrXBQ}8rnH*WJh^jc-b`n%qN7q$5RM^!hFVqDYBey))I&oa6cy3$SkiCgm7-^W1j08<TndoRbFVYP-'
        '*57c$Pjk(Ab;o)AkujoOkdta+rTW8O=J*sZAJPryNNg?jQw{=g8gBv|PMg&c9S&l5%;C4oCR{S8c&GjCH;4(f54cnNg|0Km(|%X=)lZ_~749p~9=0XA<WQnm'
        '%Bs3&4=X<|XQSLv80(Q~RI?&S|Gy&t@D+K^$cHZ<`LJ)~ug$eAvV-!QqK)rZc>`O4wSj6E2HJ(f&nWIQnWp~2zGelpVT8HV&>YXbC??M&kPW<1e23qWqE9n;'
        'T{7Rq+Geyy{T1N;ntlKJ;tZmWDHUK9z{~m@VM+T-cHW??;OiQ{e0TGbLO>HK!=|V&!*zIobnY{idquldo(4@>dzdcAPo5Y}ty*2GkZdlxOphpLU<tOC^A4Qq'
        '>}jeK?rjL`OiDVks;04vKzUH5ylH2L?S&3|{of}M;Mo@t-bJ)N51VRHX}HEcH*9T9YbbMXQF*_1xVK9jY5T+5^Ecm~9CrwuJM2G$uqleRN{{^UuO~h3YP?*C'
        'EOV_7vdHjCF0|AB@))66Wk}8@eb6o_*RcO`bE&0qb-Ob^ZVTQE?ORu}TxI*Q09c(7?BMN4s-@rEiyy;>mS_GLAhQSkmp+_YPn8x@l8K{OY4ef8qE?A-KmT5;'
        's`W`%XyOS^VfK1Aq#4;}D)+nP8GfJ*wxc=hMe?RKR<txXJvYl{IG>B^Cq*z_{v#rDH-G&4N!_A%{WclTNd-Gu21mj)E5NTw#7c2u5y6g>Ebyy!W}EKzjl3<x'
        'DF1QxKr{G4J$`ZX*xkv0fh?Gdh+9RW+ih*gg4bJMetrJ!^Sq6kQ2V<OPu`Jyskr|mHWpbB`>(mL;8<8WzIAs2vx$`ybg485L(CJ0;p3D&z>7z3>u%F1SbA>u'
        '2&ry-36GZEZyzn!=DJ@t`sL!;GEvKeqhbx9tpf2(VR*hE{9Qt@RRFd*rjLz!sGoCn2frNF`U=eP`Ps=$ee&)ciulh=_DO7mzCC67yS>LQGdVqk_wAF#p37%$'
        'L@&2h8PluS>FoT9HO{N8DkF7oNpYlc{Qp8$IyFeUgTVmE)RutU#xSIBONQRm`sBssJz!d&d`_yNtXJ{}{|nyPW_b'
    ),
    'figures': (
        'c-rMX?Q+{VlK=G-T<%qkv@6M$KQxJEb$5=f%(^DyxGZP4ro1kO5+R!#iqw#l9j(h%_Yn6O_ZIgi_axWd07w7?_2DG5b#<yrERjG1XmmHa8;!>Af48^H(!EI-'
        '@3Hv7S*-H=BtC5X-s$|I<4ltvjPC}{GM{x`;gg2zx{V~A-m@%E{X9v%S$MZhS=L>w8bkZ1(KwHDmip5?d|-}R*_pw={GK_dQ4%L-tJt51)67{cCs8={^Dv1!'
        '9J~|xE0%T}jVrcDvJk4S@KLhN8y)l0xV|@92{R{z9>f6)oIG(N>_!l#Y??<aX9<r}=Oc`i?R{)J6PASmqfX6|W!ms#KU!q~VG<?NPnqM#0li#ANuC34`13)d'
        ';lPK0&^P?lp-H~q?|IWCRVwTsc#DiJgTzbOjHN7|sxJ?{2R{k}LW}bF2!O<SnnaoM@Ys|6D<4{$+}P}zrYp@i#?z3+c{mHF&LjcE`su1VKQu#t&_$ZurT*Nx'
        '82x2*<)kcK#+j3q5G3Zi&Ux-6F-%n$KO~=6M);rlaS{WOL=7Bihu&w?G+g8`M|XZ0XSw4?Df5GsGYbJ!_z#$Uoz7rBmQgN9AxT4+9wHvX9ZRX9dB&m#hGlaA'
        'ici86W;VP-m@1_@86$*e1BcIbBZDThHjq)A!=EUOne!M1`91uYCJ_)9e3|>{Cx&1A&w!FBoQGI8TP_w!nnP>Q#rq@*vNjN6oaEFJ0%LbT&d^bR%I1KuwsUrQ'
        'PP7@YnU8}4)>tfaWU#r9C;=XV`cCq|(liVL7B?m<fDetTWC|R1CIA!!{+OhnW>NBJn@NS^Vr(gLU>2Qun5GFwAD_h>Ckpy7y$^FhH6TYYjG)JYq{ca0I6Xl2'
        'f~9!^^WseJeUuz54gW!d*NYy46adD~I~kycIeShr$yoZpMaLYtF@o9k@8TqbE@XP$z5`>9L6~4!iSbHvjT#QoU}ettARc)U2S#{Ll<Aa(8x=M(>c8&v_HnWv'
        '{S@X<5U5e6*X#5SoLRI4404jQ7Ib=t(1$q;2fC0!KVTAz0@*5tSpj4?Bu7yaO8#v-kM~Scc($AaX|5awWSM88|0S;%@KCM-7w5oRqLW9RUf)TR$ui4hmSqxB'
        'gNQ^aE<s6dLn*B>wj|}QgEcd@@UbE_)9D?x8~n>LQ9b|OFF6<KowF=7B}BX`6;-(vV7Vy!*qE@XzsyMZuTVTa;7Z^J0g(uG2n#w3KW5PH$&#~E<S$^56QDp4'
        'Tz&*{=rPRi5lBNF?lJ{Z$c1dsc!37EqAs(~Y`VnaB+qC**yja^m=Kn#M#_F&0;jT|t1qdv1ObL3fQMv(7ztxw|E!*8s|@`1{AzrS!=5j40StsgjBBYg2PQR^'
        'Jdn#L4}cDT5+Ro+;gk{D1!v#%ZXIN8NBamGX#?Y8C;PnU!L8(-N>Sh|Ugq%DDVe9rJ8i)o2(@>}s~p`HL|#m`g#kL^otPR<+c(FzC51ZMQ@aMEU&8PNkvWoN'
        '7KljxkV*PI+Up(fAr%vWaTq7Ju#80tSKvqJBE8P>4~_)PK2M`1XgH0=&%iX{B4V-_)CILhOeqW64(r}^osYftVf*+a%;F5Wub`ebaiO{>H2e&pyuWz&3+ZbM'
        'AC@o&X%Xae;LtKsLNADe)_DIR2?JReN6hEB!B^glQ+PpFMD!#BfN<JKcd7rVECzMu8;W)!CvyU_9fts{@}4h+Hhr=2Cq7}FG)eM*<*Lp8zEIXYO@RZz@VNW;'
        'Ovvt<Vsk$ZXYi{#%^s*+l-!|eWYxZamlyD`9HD2fKStR@kFR$gt#_WXwiyNcUb(n=VsYE6UCi{#+5Y$5{wuFo$GmtMMSGSt%&hj>>j6^uB4)hn^}ODJr>tE@'
        'hkJ*fT(dl}W|ajhT)bCVkjla(GkRj_Dc=cR*>kL)n^Nc&lJ>UI7~v8K{pkpmYLS3EUe!p}##^9T=1=cSDut=K!<8?9S%&{&kmOQEj&GgpQ@AMl7cu{V+a63b'
        'Dr$?SmvNpff!G3)b2#=eUa}tqw5wa__hSmniZbAc!I~);#m~)N3AhL?MU|4#5*y~gH75XnzM$$xXT-Xkp)EE3>HHEzNDSjEidRZCu10^oI=>!uou9)jLqN*~'
        '@CyUwzc8<OPcl-7OF>ZJtmt@yL-#Z1&Aa!fKaI|Q0IZ}eBB3l6S1Q8g9`rjRr4&YVy;_h24V9%6mixLjl8`9OVQ$&g8FXOxpy2@>V4r<NLIvd7<yt5{g}Q?z'
        '1R^0#g5@;C)ismg$wWfrv>P8g9WH4<@^J!-3!G*CoGH*t{XtOOHbZj$R*;ZvM_GQ`#mcA_1Jmj*%$ETgMja{Op@I*Fk*B`&<V&@-B8TN7s=Z&qqDpE?xe|Va'
        'h4G_`$Cy^YOj#a<Nfd<>qU>lTm&vU0(OjB8cFDQXn1Kf5c;0M@x|!!W;hZc3#r|qE<g@e+v{E1;@f&5n_<HYW_t32TNy~=#o@DYH*l#&nX!18yVV}z%OCa9j'
        '*Kz_B0XkMuIa_7CuZZJm<YyTSUA!wE6(vCB0kSAx({EVo;s1YtlE=&8(hi-L#V+9|f03_9UlGsG!YNps=Sa~&o^9tHh=;cGGgoTb&ie(|m)p)5#VD*VZRdK4'
        'YyoTW)#;x{<Lj%F>-Se4*!|=4_ZQC4aeLifcfX@>Ca1biF}t0E{eymI|8?i!5ay^s%<sLwIvZU<^(GODi;Jt<=HK}eCx4Z+Tf8^d6#1g=w<HsI-jeOVJ9#sD'
        'H--*8DIFL%?kOt1l8T$A(%Uln2o)~rTc=o%&B~w#w=83Untv&kVxvw~EpF8h9w3c3YQ_~5L9O$OMYmT~hYK8_y4<S8EPy++Ys9i6Tr&m$&yB|E$;JDN^V5@e'
        '9!<vh+tEQW-{b3(A4lF#qhG3X{zP9quF<(|NBMcE|LnWZ=KMd2uO4Y$#d~e0pOl{-sagnDc}6Oh`RuT4%rfUsmZ#6vNq#D)s<ePA1(j(mD)f%TQ;q`u?AY?G'
        'FP}Zr+@kyWYkmOUX|ZldhXV*0F3*L7fEGl>N$}M>ary53^)<*mop>b83`pKKU1=p?H@wXj0<6|<t8IciTm!D&h_J)zuYqN=GgK%VT^C?17Dxe7SRZBq6=;PK'
        '3c#9mqJUH_6Hj{0YCv*6+61y<&};;IxCy4k2GUSzVPGq~inT;7R37V-UFGir`jG*NkN?O6rmY%V7O15Vuz{2|S&4=t>O~4e5}iim@=`YF$Gbs4_7oqn7bJ5('
        'TnoOx8+c#{>9kj|qc?yO?%n%~AJ4BpoSk353U)Obe|UH8!IE}%es*#_QbYtng~pq1ud_xBpfFVnbSeepzz`tPVzbb$-NTM15kOnekVKAjRZfUCQme~d@Exez'
        'sW{;~)w*hL`ZT<Lrz%P;!&EE2ri!mcgDgNfN&7J%)z~yaoiixeJ`HENUY|q<Aa-9fUc84)iXT}QJs?I^F(<vwG)-WF%a(02!Kiz{>mJlKy9T5Q>X3sv>=dFy'
        'ZC{8Hqs|epbJVHw8JL1ga!MOjj(L^i4iGXZ<o+a#!hGe)T|O^kG;<H|(5TthNv+pc$-~MlDs@3|auNK4c{^59iLfZ%sTAkjcsU+@ID5Z|4ZkNOeBUvMLwA&U'
        '=m!Ck)zR=u$HJ=uhB8?oRrpo|%$Oq3y@em96F>f>4saI$+DIG!UfPEAKWN(i9UD-QPi!@CfCRYX(x0%X{K=gVo>nAT=0TWp!6-JsVG)Rf0hxq1Q0$f-OJ}hs'
        'ygUPL58hzCw$tr)d9g*xs3cGlvz0-$aQHPp0Pjb1dV$S@wl?67WnaCGMk8P|Cr%!{FiV;(r}M9fzkxsu@LF8h<6S7#lqPw1nNM3?Ka;<ktu6rI8@f~|Td?WS'
        '%@|I5=C(vb*}Z>o-1lZ+MK#f$=d+_Xdk=Ii928w%zS*J|v~z{Ey2*mYO*fsmE%deB!=y%nNbrOE0Njlg#tzI@+Kl}9B=85KWEVHyo4x&mBj*ne{%N(HiR-oo'
        'TJXSqce%i!Hwjcr#mIg4o_!9&JC^0ma%2qH>+eYzHw!wG9tNT{eUN#&Dw`8v$-&Yk5Ll!6xkW?yjKdHxZ+So1=P6s{&Pkr9;bfV!5pHQ$;JJ>=p>e>!E)y6H'
        '2V8|L$HBMQJn2=1g%!~!nT-a4Q^p&6yoZT&uG6K0ecTlp@*TeB%zdJ|&tN>ax==q{h`FQ<s~PS3Lj@i4ixfnE(><ZNMhB$i(5~Ck3Q8Yy>5In-L#_$|*jAT)'
        '2BSP!g0NNE)o`y=Bhms{oJ0?Vr*N7NacpHxf%msC>Wd_Qi|dPl`m+Bs7rQ|toMAyQJh6f2O=A#Pq0WRE;lE{v4_$ElB0yxdZ3BWju&h4k0Lc@6US-PL()Y`G'
        'we8GdjpFn%a3)C-0p{NNQN~_M8_vy>XyN(RN+)7^n)Yc$gdto8U65UE71O9uTiJ8+Ic<sNGl}5c=du>nE0%e2V#jg?r}%((EnCiJ;b(M7-u$u`-?}&`8B;^!'
        'cW+jE7q{+?6Y#|0YFyfll1Ghe>bO;m5j)wG=vOOYlOSfE+!a~@D|FO>)<OS;^#6bVm-}o5iCwosB6}iEYN`?98ip{rjF_Kcz)uFFa#px~8+VIk$s<R(N0Qzv'
        'yt+thrPSjKK1~&__d<-b&pDDNb`3zfCp8F!-cwQPzK`Aaqe%QIHx0<=RQ!0Z2=pAg&@^3(wicr3kA+8=8u65$?181L{GL>Aw9ERnO(~A-7T=$Vm-7k3jpPMt'
        'pSYNW&L5p#Yk=w*Eh=G~zCu436sDaG4bC&c0sYU<B$!fjls;<3NoqKwHJ7RC#VMRI_ENfIx}tw-P7H!j<9e-kQLn9`7PP%ryv!D#eg%xTxje^no_qr12atS|'
        'OTy%FtvEqQFkLj{{3&a?KOi-<k$}On1D1iWL{+lhZ9@?Q9Ih%@LSV1I_|cLpBTZPH<W`rAz94=hF($(?pj#X&bqPUIw+z@wPD6$Ax_1n8Of1wc>8O0z(!>wZ'
        'V&+_cj%r8?v&S>{L;MMJ+}PX$1!+%!f0wEPU%*CWgmEz>y^Uhz{sOm(n@~YhBsrJKQicj<N~3xyBxKy*Vi$PP9lD^_zIxSCUyTWh+=3SKk|(%l>uLpUmDj~f'
        ';m#2cb`*d6u32||09B><##kyP-mooW3M){o;Z<8~BAIZkkyhvSNGpouwZvOhfLcSZ&enxnIZZV(YtKo<E@rC&+ZvMfrrdLGnju*Rh`+@}gjPX%JYo^R!ipv%'
        'v?xHU;OI@dh9E$kQRtyn%d2SRij=5dpDVg~$<SyT!b6)<2ah*d5-Q&STRX&B?ZTpxEZeFWapH|M3&l^5!amSJk#U1Evzmu2FMJ&+wBfhFkY>{g0IGn@q_t}9'
        '1>R8MTTSd>dC@<h*+7ImkPeI?4pPwv8VZmPnYFKSBb|I8LJ(muV~p@XyRxipYVs6Eu}39Uak@qqR}^VOwHI$tzKbEDwnJfnT-7H(1B^gM1n~YEinR<K#E$dw'
        'R(p&G9M^PC*E^N*khh~nCI>x2Uvv{0z{4PX=TAdte*-{OAr^E+<Rqh-n)GwPy=BJGL|qdrRHJM9zRD#?advrhP@&oeUyHCM9%aEXD@b83f~jDMUvTu&GD)p!'
        '&<g|?*u7Qz-NbhSSgS3+iEd8slLX{K;jAF=)s7Zk7ts{&1@{jYt&&NA7IiF!Dhr?oiAb_l^|0XjDNomk0JEWWZd%(d)VyIb%Qh{`Ggrp`iSV=m>f>d<pIyB)'
        'pd(+Fix`yt?)3f7m+wZ`qwf$gWk4MolHZ+Q6mOr+6+(sv5}^P2kN+$}1H(*-e>I{4bs$PSAP1ob16lx8!T<$8Y5#&zJa(%Caj&+uhEFL>yv`f3y**4oWxSlF'
        'R$<0+BN*^lW61J%)~YF175JmCK!-IzWGQoL8{zhfGe)w%DlK6EL_2IbV>N!vHU<rQrSlfniN`>GW%0v6##(VLwOB8p(SeLPrs`c-)LEQ&=bw<f_>YVnZfzda'
        '=_Q{O|64(Q|3YhV)!*9E>HP-UxpkgoYm$h`iClDnJCYw=>^VP9ICSYaAv~EZXJAp{qaL}Lp^~za@o+cf_9FJJdr9E5S%^^9$_}=5o#QbWq5%zmjW>YBVGvG1'
        'Pl+FlS6LnaKMd7I&o4)nGAvEiGUMyB_aCnHQfT~fY-^8U)VL}Qd*$PjdNGU-BLf9k-WfEnkP~PLR(<evs0<fWt+-7e!;&}!!O%5aYJl`E(Jz&$Xj@oWf`9bE'
        'YHdkErb<P?dsC&?G<nR#CJ+V-aym^We+NbFmeXo!EzUv~1!O;$MmR}X)hZ`>tv<*|)*4SzR0I87eG1xF2!w1t=7YnC9}NWUaj8$ZWJG;Q>y038-ETk?s6Ltl'
        '&dPC!J+9r^(%h6mUW(Zx&K68&m9MNjHZ)l&V0oH7ba77pN)MX6T$|4Vs*Bl>Pe{wC%4y^t^b}3nA{D}HBN5W1Yow*tgtyJ7PXjDPj)YGjQ(f*UA&r~2`rP^^'
        'yL^`-<ZG`CI450{6x}7A$NM(henOM><iWJ2TywqB3?P0e*+fzHk(bRk=*@dp^}m6l8S3sHFhYP88AiN2GWHC0t+}RO%dHu;V}K~ASDs{E<~IQBR_x%^ZKQB5'
        'Tqgn1N$H-!qJ=kA4#046RrY>#dGE?Y1XDetj}9pC6TkzZ{=pw^OwPb=?>JA3XSY<8rN(OrVIlVLsE|H|l5$WOWKE3R@KjSy`bg`hRX67B6Ks?zr@ENU(&WP-'
        'tdAk8PPH<Rrd3*uM{T3EuJO8<iJh8NuvTqYn~JY56v%o(zXT6@`F28JvJ4~Pp2#>Y|DmlM#o)VB{DBY?hK5Jvz)K+U8E*}$X$#yHr{=sOMng1WdhJv?7E3Z6'
        'Y%?E=0~WHLbgc3Hw046Co@5bWR**O%fZVdS-EI6lpgRfmc^DaA8KesrKfiZA;HDDDnYdG)L&6*ZF$gv~tmRMex05BBO^p~xQYxM@4D!n;;5X`;3Jz}{#F6+?'
        's#@H<;Zc8ko_7lD_Wtabi<6(vPsiTR=i~AD#gCo`i`X9Whj6^pEDuyXaG~7bqY)vrP=)(dr5GUID6kH%Ww5P_-#3Q=yGrX)*eb`u(P-DCD0FzJhpg*LA?jO0'
        '7uQYhmKG=T(T-emI^Td&#tf}%iyPBJ_+x1AT%#J^G|ZEFoNwBd=Wb`|R-d8;AQ|gM>&I=&2BRY-eMeXORCqyC1ZKAEAcS9f4Fnjjy{f{wYpNUkZGgIM8-}Xx'
        'gBms~Fp&iJP+h(|@ut&Y<}HAyS<qbsGfNN3(r-j`C`G-c7<Y+*fUAx5PQ5!HrdJ7a1gQdZ6ZjPDy`lD6R`BNJUPrIK1(1O$%~58-OF7Zlw6#?pIuClarOUO|'
        'QF02@Ww7qmR5xb1b|hS0wrYgT;SNHk{&2t-N|?i~B+MaqcI}YVV=Iw!xU0w!5{65ZS_yK95~M0IR0(m25~3mz3JLIHM#UA)72I%e=J^v0#0A!&C}O}2&}EW('
        'ly7D|H)()II-rs{#opKspx%10f|jQaFo+e#6|hBq?@?$r9dB3%j?T1FdK*7ox$?=l_tyJ$>dhHOk0J5W;a@rk^Obo0kTdHPa|PwA1(O_6qV?A)n()6aSzgcb'
        '+dz?{kY^{VLSSni@9)6lwnIB#$ld*|xVul>y+zuct(dz1ZJ4^hnW_6<#?(9}p@!^d?h`aIe<pISltEhzU!wt5g+LmF<rKw<e*Suj+>9G2la&|PNVpEp8~8>}'
        'yz+!j0qdUm`P^UZ0PWOc7~G7$aQSH04xN=TN!s}UaU8l9Uh6A}i9ToW{eya>tXB#ofTI$LLJ;tNZ&$m&4h5u^$+I0h0DV<6HM>cu>M=hdt;)PUc7au`SO_m`'
        'a*`JbFYRbzz4Tg`N3+X>UnI<mY+wA0P_5ir)7%#cE)A}bTy$8aDz?6tVvA0sRNudn?|Fy5mYO2+O3)O(rLTu4#o?Jn8uDX<n-mKRl#w8l$iMTlC|O|Yimerj'
        'D!lR0sjzy)0W4rkBzpkA(I7$c&5~sbOo-=e){a#zxSpYQkVfB^!ep6(z6G)}NaoTaV(|k*1g+=(SA|3aBaY1{5%e)%qnuc(9(j{6j`1Yk7qs73&cR|W1K5kS'
        'CH07GA%!s7ZH*XG>nH+9Jxu%}$x}VZUn_l-PAiE5QSg|_0@sM14efoA_^H5FM9}N1h<UAwpx%D906_SSgwW9mM%_emCSM|fjz9u2VC6Uufmb-Kaz{v*Sq+kH'
        'gwYsY0Qm9o85cvS8+a3-NleqT!2;Xgr=2G5p#_E7^FigyznYMFEf{p)!ypDFBib6z$NrXecI91)PTiS#t9i)2l!HaujqN#@XJjE3wxsghj*V@Nel0fUuoX6j'
        '4p-T<r!n$@rNqH+#mL8kk%?Qxnq9atUc%PLBZ_*-`FY+T%=j6=@QYcxe<l@M7KL}AJVd>Pit&vZyVrXGXCHhM&X(yrci?JE;*2lm>cj0>8WZkvQ$SC;yCo|h'
        'ek(pc)R_3t!lH-26^r&y3gx%p(uZs_+1EfDp_cg=X6(z<e#w^KpA{y2K}V!TXa8KLJbVFD$~woU4*e3X`O?f0EsjWXCk0Rp28395lR4h0oVyD<S5J7o$h2Q`'
        'PRwNHH5swjS@RoO_%f5e0$egO$Ed=|%3OSHPUgI(CZaegX=;%s811uxL8>6$cwszEIWuA+Nb-Dbu>QwnE}wz5>}}U4yOFsKw3_aHF-NN`!#P@MbiILXR_$#)'
        'V6*0db#Oa2e^^$Wh55^D@iKuDfVAdF$5UQ1aeTe>DB+Oos3xO&;ig747goLauw}GYPpnbEL?s*k+Wn19&HvW%64n0H#7*$zLR&;l9My$Q@LWS%Ksa6tp}z%$'
        '{#po!TR=FhgCJ$_mT@0nD4}1Ss2dM5x7%j1{OZh55oo)lNR_lfCaXFXhxdBt8omA9tMs;3$u?Uf-&z+tdO@Az7u4y0p&=H|YD4>s)af4P9MWU;@l!>cB>`x)'
        '3<*$fK^>KMAyjg`R_Y1l*1oCc544o6**#S(jy+Ct-=BD3&oRt%R_00Dp#SjBSnQ5HX>$Pgp}Z#Fc(9#0(^^Z&Yv$m1yE!;6=Rl3P;VBX~O~-Nl8%x5)FPRUJ'
        'cVt6itome`1}uK@B)tA{YmN5B$nS55{C<giTeiAQ=<nCR*_14A=LvXyQr%pB@f7s-|Lk2~e}!qd)sCc;8daOCYwhCX2Y(0)G@f?B^AW##;*5kU;|8Y6wORj^'
        '$wA<|fUd)!Fl___D#DZ6n<7ETcQ|9;wH^Mq^-M~L2ITk7=(9h~ozr*c5+EJgq4WW&1_qr_)v59)=|<fxduHiDjW)eTt5&sZWnQdhxN6O3%vwNEl{dkv)020x'
        '#mhF(1m}p7J~S0kb$&D@N~n`zxEN1MbZhTXO{5UlZIDRmv=t{-LANbCBtxw+I$NhLv^GY{G}MC`oD{yv*$kxbJqKmaH#(ES(_Q%KUDkyE9tPzp3hgK^2Jhoh'
        '06MoIuY2fTAj>U)V?}ZpC{(;VN$-{@>n`z^a&#$RoR*-Gp?k$B<9O}0mqp@yx_zO%ztFrcrK^ZYo*5FTmP_YhrAOjszwFN}TD+BoMQ%_&j;<dBh!ZN_bnSO~'
        '01ab2<B%VF5~r|-grpiMMF`Cc>sKP;ckY0ZizcIcMs|mzTvD=Q<%VJ&%F}8BPMMyifTs{m3AcQD4-+_4lIG#4fK%|Oi++P%aV`#M^%KUUIcG8jo|)`kQqsoW'
        '^UNsjR(&)M?nxFaFSIeQ!r|w*JMyx>fG=))JTX<)>N@&b<etx*WFF>ptC{6m9Y1>XSDAX@PqzfnT2K1IB9bK2DYLuzQ_jR$jxnQoy2ShdG~R;c#vd*(uSVl>'
        '9q>5$HgJP<)xkssR$z6Nc*4N?*4;5CqTo)sN6;e~HsR31fkFgo{b0YZDihZ);hl<^KVwe1M9t3CPWLUkhA5WD%Tf4mFPn?p2MCAGshU8ArHhZcHRdj5z@P$g'
        'Qf~@Z6S3hDHA<^N*Co!JTI#`65Z2jml@|ze{Zgb@_2^Ko=PBtn#p&_#ykK!?cc|#_Gz#stmeO;55>|;O9lR|G0iBR*mfhxCPqLo7rk$2{;1N?6tvxDTc#q#y'
        '*<8KQYqBcRmEdIFlb)J6Z$@w5Uyb-#XP*2m>uMAs61G+H*oo)E+LUn(i3UrZnxtB&94d{r%CUi%##p2d1rg<R`PV(&*6Yk7|Bhd2MuY}DyqB_hBCd;*$;yh`'
        'APF)bA#WYCWI&H7S%zqKYmTVqKa@H7O4!Ee>6tVVs`d4@tISr{?B*spQfiwlriN?QWks2Mld&M9XnEA>;#J8|lQJzBPZdTX_#%y_p6*i3Kk3x<q@}bclJ^aZ'
        'teRz8N-6as%%Mo&X)Rc6WDac`s!A0#jeu7isvXuHdHrzlVLUp+6R(9dE%8;(Msbu}6G(+nvRwtYgKSa>$9k)rO-vP{%9eYZXOYg6BCS#2sDKNldr!7NdYhRy'
        'Ko<cj4k`%V!we$Tp1THCIV~uO8p^ZscFYz*^Iw8f9mrR~W+jAv2Sn<FwxLq*u!%&4Y1G_e+ObMIHjP3rM5o!*%kf!UI<FC@au{3OwYUMH)v?=<s+6jbrI09t'
        's=8ng1yP1<WG0g3#q6#DV#mIcnlh>pMp`>%BD0b|>w`|}x@pPfZgK6J=}$G9uliZcK2}SVsnv&+w$Si|YAV@an*nTO+qFZmGON*$m0hdI)|3I&$!et>wZO1a'
        'XC0yFon1pPVh+$(E!!ywMyTKAxetX}6e=`@cMVuYylo4(!hUt2)yY^7JlhNn$)h`c44N$~41H2Iafo|{$ID{w>f%?vf<~p^GFjIb{cSx{NF7D!wNkMhu2U<8'
        'MueG26oaHD@sF1Br={d`q5PC7ze^!o;X2Hw%N!G<Dg@b5LQ|kZHAPGbZ)k>ZNDFBwzbMWa=;9FnYu5?<VQo$UyW|~ebkw%Z_)z&?7O1f(hQ^GjUYJKe))V<s'
        'g5J9Hz6J_!xa)zRkmXly%a@_Z=Gql$^DFI$^emgReAjwxg>TSNFI09nCm5C4kELZ%Bo1vm-|^H{Vyw@l&192W65cyqH{kMmfH@ik)mGqCv`l%m3b#oljG|`>'
        '(N|zul`RSI%70X8u>fe$W3_`7fmd|V4A_)iR4LX>(t0U<qUf|xQnhN{2D1|1_~fEW&F`uZzFFl|ESg33aXGV+pubb&E{c@llo^&X=*nV7ZrJfO#MN{bPRrc@'
        'fO|3e%jn9%*u42-nVY#_3)(Z&J8heTuxadlM#&u*vfI=ael6T((l4UWIAvhEyJ&OPUf6p8k;iYDuw5su)RL)eDf{D~yn{<k`6$JU@M{k=7m1eh(PFVtPOEZB'
        'XzA?9XH)l!o8}6Kpspk)Yg(ZgOe4UUNr9}KY44R)_RHR>{M4#)rX8ZnPZE&k(bKYY7Stt*Do@mcHH(VNbC5L+c2a8=tGL1=1IJJrT)xzL<&q6$%;i@WD#;L6'
        'dGSkJ?a{i~ztxbzFM3C4;%ira(a2G(=pMh74X;wJZoYT5tN7IG-ilpoqu8x(n_8C0;{CIFwFzfiPh^NR!yooVEQhBtr@EnR0^ylJJ<N;S!~)%5IbUQ=UZM@X'
        '2|$hv4=|-^hJgKS8iqr0O%u=OXi0Z!tlcA;JJc)z!|}i^lNg4<=XN>p^6ngt0(zC^!=9Fe+z(c5CDVFs)^<JrG(I1gGe?=rU77RSEwmCdkGagM{o!KoYV_Bu'
        'b11t;S1$n5wFILxx3Ww+&J&D$uzwIdcVz-|d1u@cRW~kOspRF92vqZC-VBa!2URn+Xj98jikUZ!Etnu(VE}@;s$=f+v#Vbm82znt$|E$gfvYE+HHg=KsTcry'
        'Dzv`m+rPzt2PLI7e<kxqs~8X)Q~M1yp*{6^#VofLXxELoSIgM(A3*q^@K6$yazzwYpnr4A?P6LX?KksHk~FVg$ESCV#Sbv?30N5MLzn_PO_Rk++Pd+>jeGga'
        '^`GBge7Jb?;qBYe70?E5nYleqA~JFFZtN{<QN+IvXKm7FB?;nh$BHsfb9Hp2zwm53k~chsTDtIF_$8F7wy?EG@r-DdX=#@W7&!Ci<>XWP`%GML@<VaK2`-4O'
        'HTh@_Vix>;tFthXwdjR;S)MjCwGmrSHBoz%f_qRJPBqbeJTtc}w<dB+rtAhhkah(+m-G*D1iAI~(seqVR#?sVy}7iK1{}Io+X0>yd|-rzRTkNr>PilA4nn;{'
        'qOMYrVqB`2{AXQpnSFfV)QY&8nlefGK+<xx+;gKoBFsnvuL+g5PJuP|0hL?UQxS%^=FSm)k6Kh1ihCUF>%h(Ur}ImLYyPh=3Jz-l#4J(q&t{SR|KA)^-Unk|'
        'odSyyzk(hE5YYLIL>x|8*Zrrllg=RWjuP5gx8BZUsFnA4BX&2P!&=7l`Dp0Yhyb(wod#8nuj+0qUsvwysdK4`9;N&~vx2Nk1!9UmUFGr>ZSu=ISK4}14c)k+'
        'sGov)v?es&e~S3tFY#=u_f_?MDcsUx<HfDhFkHKQd;e3tbfb*}(WSbEXLRg9#2kB|>rFq^c0pdG08@&1&X!ZW7-Y7LqLnjQ=0yav5`-unqaq{>1J&iF05t=8'
        'dvgBnwqmJRlpMP?q>id?ATx6$bjgNLugNp`lun;_=wR0wGp65#)7XA_d1}05S=nw$DVG!e#a^OZKEw}8F@52Xzbpl{H6xwcT72<>vGZcF!$R@0MPic*^=!H3'
        'XD2zYB*T9TRNa~RS<b`2(0K5{;LP#-4DH!!24~liMoAGMCan_%<|OrFwD5>E`F$_}ULhk44w7U+BCd5y7S1xJl)F$!u=ZV=0NAc71qC1ePLpfa+p|D8aS`=9'
        'l$zxFIC$A^Ol1d#2hwvY?OHwh*=8p``mM2);>I#d{Vb*ng_Jsx&ucb-N=v>0!}mP&cY7W!)t*OZt~`M+9&S2bVcJ*FKaf~}C|j+@{{lcnFq8'
    ),
}


# =============================================================================
# Embedded engine loading
# =============================================================================

def _engine_source(name: str) -> str:
    if name not in EMBEDDED_ENGINE_PAYLOADS:
        raise KeyError(f"Unknown embedded engine: {name}")
    packed = base64.b85decode(EMBEDDED_ENGINE_PAYLOADS[name].encode("ascii"))
    raw = zlib.decompress(packed)
    observed = hashlib.sha256(raw).hexdigest()
    expected = EMBEDDED_ENGINE_SHA256[name]
    if observed != expected:
        raise RuntimeError(
            f"Embedded engine integrity failure for {name}: "
            f"expected {expected}, observed {observed}"
        )
    return raw.decode("utf-8")


def _load_engine(name: str) -> Dict[str, object]:
    source = _engine_source(name)
    module_name = f"clonodynamics_embedded_{name}"
    module = types.ModuleType(module_name)
    # Deliberately expose the MASTER path as __file__. The validated child
    # engines resolve code/01_core ... code/05_plotting from its parent,
    # which is exactly the desired single-file repository layout.
    module.__file__ = str(Path(__file__).resolve())
    module.__package__ = None
    # Dataclasses and some introspection helpers expect the defining module to
    # be present in sys.modules while class bodies are executed.
    sys.modules[module_name] = module
    try:
        exec(compile(source, f"<embedded:{name}>", "exec"), module.__dict__, module.__dict__)
    except Exception:
        sys.modules.pop(module_name, None)
        raise
    namespace = module.__dict__
    if "main" not in namespace or not callable(namespace["main"]):
        raise RuntimeError(f"Embedded engine {name} does not expose callable main().")
    return namespace


@contextmanager
def _isolated_argv(program_name: str):
    old = list(sys.argv)
    try:
        # Prevent arguments intended for the master from leaking into a child
        # argparse parser. Child interactive behavior remains unchanged.
        sys.argv = [program_name]
        yield
    finally:
        sys.argv = old


def run_embedded_engine(name: str) -> int:
    namespace = _load_engine(name)
    main_func = namespace["main"]
    with _isolated_argv(f"orchestrator_clonodynamics.py::{name}"):
        try:
            result = main_func()
        except SystemExit as exc:
            code = exc.code
            return int(code) if isinstance(code, int) else (0 if code is None else 1)
    return int(result or 0)


# =============================================================================
# Console UI
# =============================================================================

def _rule(char: str = "=", width: int = 88) -> str:
    return char * width


def _screen(title: str, body: str) -> None:
    print("\n" + _rule())
    print(title)
    print(_rule())
    print(body.rstrip())
    print(_rule())


def _yes_no(prompt: str, default: bool = True) -> bool:
    suffix = " [Y/n]: " if default else " [y/N]: "
    while True:
        raw = input(prompt + suffix).strip().lower()
        if not raw:
            return default
        if raw in {"y", "yes", "s", "si", "sì"}:
            return True
        if raw in {"n", "no"}:
            return False
        print("Please answer y/n.")


def show_intro() -> None:
    _screen(
        "ClonoDynamics — unified analysis and figure workflow",
        """
ClonoDynamics analyzes longitudinal TCR repertoires using paired technical
replicates and a noise-aware, replicate-resolved architecture.

The project separates five layers:

  01_core              Genuine longitudinal preprocessing/state/transition data
  02_pseudo_reference  Technical pseudo-longitudinal reference/null behavior
  03_validation        Forward, fluctuation and temporal longitudinal dynamics
                       plus support-conditioned synthetic positive controls
  04_controls          Structural and operational-observation robustness
  05_plotting          Read-only generation of publication figures

This master file contains all orchestration logic. The scientific Python scripts
remain modular inside 01_core ... 05_plotting.
""",
    )


def choose_mode() -> str:
    print("What would you like to do?\n")
    print("  [1] Analysis")
    print("  [2] Figures")
    print("  [3] Self-check")
    print("  [Q] Quit")
    while True:
        raw = input("Selection [1]: ").strip().lower()
        if raw in {"", "1", "analysis", "a"}:
            return "analysis"
        if raw in {"2", "figures", "figure", "f"}:
            return "figures"
        if raw in {"3", "self-check", "selfcheck", "check"}:
            return "self_check"
        if raw in {"q", "quit", "exit"}:
            return "quit"
        print("Invalid selection.")


def show_analysis_overview() -> None:
    _screen(
        "ClonoDynamics — analysis blocks",
        """
[1] CORE — genuine longitudinal Steps 1-6
    Needs:
      • corrected longitudinal repertoire dataset
      • destination for core results
    Produces:
      • repertoire characterization
      • replicate-resolved latent/observation representation
      • diagnostics
      • multirepresentation trajectories
      • longitudinal transition table
      • transition/support characterization

[2] PSEUDO REFERENCE — design + Steps 7-8
    Needs:
      • 12 pseudo technical measurements
      • destination for pseudo results
    Produces:
      • pseudo design and 66-pair bank
      • 2000 compact randomization configurations
      • forward technical-null reference
      • fluctuation technical-null reference

[3] VALIDATION — Steps 9-13
    Needs:
      • completed core results
      • completed pseudo-reference results
    Produces:
      • genuine forward dynamics and pseudo comparison
      • genuine cross-replicate fluctuations and pseudo comparison
      • temporal scaling

[4] POSITIVE CONTROLS — synthetic recovery validation
    Needs:
      • corrected longitudinal dataset
      • completed core results
      • completed validation results (especially Steps 11 and 13)
    Produces:
      • frozen empirical calibration
      • oracle dose calibration
      • R0p00 / R0p25 / R0p50 / R1p00 synthetic scenarios
      • recovered Step 2 -> 4 -> 5 -> 11 -> 13 pipeline
      • verification, recovery summary and positive-control figures

[5] CONTROLS — Steps 14-16
    Needs:
      • finalized longitudinal/validation results through Step 13
    Produces:
      • interval-position / subject-composition controls
      • fixed-alpha operational observation-domain sensitivity
      • cross-alpha observation-threshold robustness

[6] ALL ANALYSIS BLOCKS
    Runs 1 -> 2 -> 3 -> 4 -> 5 in recommended order.

Each selected block uses the exact validated embedded orchestration engine.
""",
    )


def choose_analysis_blocks() -> List[str]:
    print("Select one or more analysis blocks (example: 1,3,5).")
    print("Use 6 for all analysis blocks; Q returns to the main menu.")
    mapping = {
        "1": "core",
        "2": "pseudo",
        "3": "validation",
        "4": "positive_controls",
        "5": "controls",
    }
    while True:
        raw = input("Selection [6]: ").strip().lower()
        if not raw:
            raw = "6"
        if raw in {"q", "back", "b"}:
            return []
        if raw in {"6", "all"}:
            return ["core", "pseudo", "validation", "positive_controls", "controls"]
        tokens = [x.strip() for x in raw.split(",") if x.strip()]
        if not tokens or any(x not in mapping for x in tokens):
            print("Invalid selection. Use 1-6 or comma-separated values such as 1,3,5.")
            continue
        out: List[str] = []
        for token in tokens:
            name = mapping[token]
            if name not in out:
                out.append(name)
        return out


def show_figures_overview() -> None:
    _screen(
        "ClonoDynamics — publication figures",
        """
The embedded figure orchestrator generates publication panels without rerunning
scientific analyses.

Its own interactive menu lets you select:

  core               Step 1, Step 3 and Step 6 figures
  pseudo             Steps 7-8 technical-reference figures
  validation         Steps 9-10, 11-12 and Step 13 figures
  positive_controls  Synthetic positive-control figures
  controls           Integrated final Steps 14-16 robustness figures
  controls_details   Optional detailed Step-14/15/16 panels

Only the paths required by the selected figure blocks are requested.
""",
    )


def self_check(verbose: bool = True) -> int:
    code_root = Path(__file__).resolve().parent
    required_dirs = [
        code_root / "01_core",
        code_root / "02_pseudo_reference",
        code_root / "03_validation",
        code_root / "04_controls",
        code_root / "05_plotting",
    ]
    problems: List[str] = []
    for folder in required_dirs:
        if not folder.is_dir():
            problems.append(f"missing directory: {folder}")

    for name in EMBEDDED_ENGINE_PAYLOADS:
        try:
            src = _engine_source(name)
            compile(src, f"<embedded:{name}>", "exec")
            _load_engine(name)
        except Exception as exc:
            problems.append(f"embedded engine {name}: {type(exc).__name__}: {exc}")

    if verbose:
        print("\n" + _rule())
        print("ClonoDynamics master self-check")
        print(_rule())
        for name in ["core", "pseudo", "validation", "positive_controls", "controls", "figures"]:
            print(
                f"{name:18s} "
                f"version={EMBEDDED_ENGINE_VERSIONS[name]}  "
                f"sha256={EMBEDDED_ENGINE_SHA256[name][:12]}..."
            )
        if problems:
            print("\nFAILED")
            for item in problems:
                print("  - " + item)
        else:
            print("\nPASS — embedded engines compile and expected code folders are present.")
        print(_rule())
    return 2 if problems else 0


# =============================================================================
# Optional CLI
# =============================================================================

def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Unified single-file launcher for all ClonoDynamics analyses and figures."
    )
    parser.add_argument(
        "--self-check",
        action="store_true",
        help="Verify embedded-engine integrity and repository layout, then exit.",
    )
    return parser.parse_args(argv)


# =============================================================================
# Main
# =============================================================================

def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.self_check:
        return self_check(verbose=True)

    show_intro()

    while True:
        mode = choose_mode()
        if mode == "quit":
            print("Exiting ClonoDynamics.")
            return 0

        if mode == "self_check":
            self_check(verbose=True)
            continue

        if mode == "analysis":
            show_analysis_overview()
            blocks = choose_analysis_blocks()
            if not blocks:
                continue

            print("\nSelected analysis blocks: " + ", ".join(blocks))
            if not _yes_no("Proceed", default=True):
                continue

            for block in blocks:
                print("\n" + _rule("-"))
                print(f"ClonoDynamics analysis block: {block}")
                print(
                    f"Embedded engine: {EMBEDDED_ENGINE_VERSIONS[block]} "
                    f"({EMBEDDED_ENGINE_SHA256[block][:12]}...)"
                )
                print(_rule("-"))
                rc = run_embedded_engine(block)
                if rc != 0:
                    print("\n" + _rule("!"))
                    print(f"Block '{block}' stopped with exit code {rc}.")
                    print("No later selected block was started.")
                    print(_rule("!"))
                    return rc

            print("\nAll selected analysis blocks completed successfully.")
            return 0

        if mode == "figures":
            show_figures_overview()
            if not _yes_no("Open the interactive figure workflow", default=True):
                continue
            rc = run_embedded_engine("figures")
            if rc != 0:
                print(f"\nFigure workflow stopped with exit code {rc}.")
                return rc
            print("\nFigure workflow completed successfully.")
            return 0


if __name__ == "__main__":
    raise SystemExit(main())
