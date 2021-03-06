function out = dynu(sys, x, u)
%DYNU
%    OUT = DYNU(SYS, X, U)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    30-Sep-2020 05:28:48

m1 = sys.m(1);
m2 = sys.m(2);
m3 = sys.m(3);
m4 = sys.m(4);

l1 = sys.l(1);
l2 = sys.l(2);
l3 = sys.l(3);
l4 = sys.l(4);

th2 = x(2,:);
th3 = x(3,:);
th4 = x(4,:);

t2 = cos(th2);
t3 = cos(th3);
t4 = cos(th4);
t5 = sin(th2);
t6 = sin(th3);
t7 = sin(th4);
t8 = th2+th3;
t9 = th3+th4;
t10 = l1.^2;
t11 = l2.^2;
t12 = l3.^2;
t13 = l4.^2;
t14 = m2.^2;
t15 = m3.^2;
t16 = m3.^3;
t17 = m4.^2;
t18 = m4.^3;
t19 = th2.*2.0;
t20 = th3.*2.0;
t21 = th4.*2.0;
t29 = 1.0./l1;
t31 = 1.0./l2;
t33 = 1.0./l3;
t35 = 1.0./l4;
t36 = -th2;
t37 = -th3;
t38 = -th4;
t55 = l2.*m2.*m3.*1.6e+1;
t56 = l2.*m2.*m4.*3.0e+1;
t66 = l1.*l3.*m1.*m3.*3.2e+1;
t67 = l1.*l3.*m1.*m4.*6.0e+1;
t77 = l2.*m3.*m4.*7.5e+1;
t78 = m1.*m2.*m3.*6.4e+1;
t79 = m1.*m2.*m4.*1.2e+2;
t80 = m1.*m2.*m3.*1.28e+2;
t81 = m1.*m2.*m4.*2.4e+2;
t82 = m1.*m3.*m4.*3.0e+2;
t87 = l1.*l3.*m2.*m3.*9.6e+1;
t88 = l2.*l4.*m1.*m2.*9.6e+1;
t89 = l1.*l3.*m3.*m4.*1.5e+2;
t90 = l1.*l3.*m2.*m4.*1.8e+2;
t91 = l2.*l4.*m1.*m4.*1.8e+2;
t92 = l2.*l4.*m1.*m3.*2.88e+2;
t93 = l2.*l4.*m2.*m4.*4.5e+2;
t123 = m1.*m3.*m4.*6.0e+2;
t124 = m2.*m3.*m4.*7.5e+2;
t125 = m2.*m3.*m4.*1.5e+3;
t133 = l2.*l4.*m3.*m4.*5.4e+2;
t134 = l2.*l4.*m2.*m3.*7.2e+2;
t22 = cos(t19);
t23 = cos(t20);
t24 = cos(t21);
t25 = sin(t21);
t26 = cos(t8);
t27 = cos(t9);
t28 = t8+th4;
t30 = 1.0./t10;
t32 = 1.0./t11;
t34 = 1.0./t12;
t39 = -t21;
t41 = t8+th3;
t42 = t8+th2;
t43 = t21+th2;
t44 = t19+th4;
t45 = t9+th4;
t46 = t9+th3;
t53 = t8+t21;
t57 = t37+th2;
t59 = t38+th3;
t60 = t8.*2.0;
t61 = t19+t21;
t62 = t9.*2.0;
t68 = t8+t38;
t70 = t9+t36;
t71 = t16.*7.2e+1;
t72 = t16.*1.44e+2;
t75 = l2.*t17.*1.8e+1;
t76 = l2.*t15.*3.0e+1;
t83 = l1.*l3.*t17.*3.6e+1;
t84 = l1.*l3.*t15.*6.0e+1;
t85 = t19+t38;
t100 = l3.*m2.*m3.*t2.*1.6e+1;
t101 = l1.*m2.*m3.*t2.*2.4e+1;
t102 = l3.*m2.*m4.*t2.*3.0e+1;
t103 = l1.*m2.*m4.*t2.*4.5e+1;
t104 = l3.*m3.*m4.*t2.*5.0e+1;
t108 = m1.*t17.*7.2e+1;
t109 = m3.*t17.*7.2e+1;
t110 = m1.*t15.*1.2e+2;
t111 = m3.*t14.*1.2e+2;
t112 = m1.*t17.*1.44e+2;
t113 = m3.*t17.*1.44e+2;
t114 = m2.*t17.*1.8e+2;
t115 = m4.*t14.*2.25e+2;
t116 = m4.*t15.*2.25e+2;
t117 = m1.*t15.*2.4e+2;
t118 = m3.*t14.*2.4e+2;
t119 = m2.*t15.*3.0e+2;
t120 = m2.*t17.*3.6e+2;
t121 = m4.*t14.*4.5e+2;
t122 = m4.*t15.*4.5e+2;
t126 = l2.*l4.*t17.*1.08e+2;
t127 = l2.*l4.*t14.*1.8e+2;
t128 = l2.*l4.*t15.*4.32e+2;
t129 = -t88;
t130 = -t91;
t131 = -t92;
t132 = -t93;
t139 = l3.*t2.*t17.*1.2e+1;
t140 = l1.*t2.*t17.*1.8e+1;
t141 = l3.*t2.*t15.*2.0e+1;
t142 = l1.*t2.*t15.*3.0e+1;
t148 = l1.*m3.*m4.*t2.*7.5e+1;
t149 = m2.*t15.*6.0e+2;
t153 = l1.*l2.*m1.*m3.*t3.*4.8e+1;
t154 = l2.*l3.*m2.*m3.*t2.*4.8e+1;
t155 = l1.*l2.*m1.*m4.*t3.*6.0e+1;
t156 = t11.*t17.*1.8e+1;
t160 = -t133;
t161 = -t134;
t166 = l1.*l2.*t3.*t17.*3.6e+1;
t167 = l2.*l3.*t2.*t17.*3.6e+1;
t168 = l2.*l3.*t2.*t15.*6.0e+1;
t174 = l2.*l3.*m2.*m4.*t2.*9.0e+1;
t175 = l1.*l2.*m2.*m3.*t3.*1.08e+2;
t176 = l1.*l2.*m2.*m4.*t3.*1.35e+2;
t177 = l1.*l2.*m3.*m4.*t3.*1.35e+2;
t178 = l2.*l3.*m1.*m2.*t4.*1.44e+2;
t179 = l3.*l4.*m1.*m3.*t3.*1.44e+2;
t180 = l2.*l3.*m3.*m4.*t2.*1.5e+2;
t181 = l3.*l4.*m1.*m4.*t3.*1.8e+2;
t182 = l2.*l3.*m1.*m4.*t4.*2.16e+2;
t183 = l2.*l3.*m1.*m3.*t4.*3.24e+2;
t184 = l3.*l4.*m2.*m3.*t3.*3.24e+2;
t185 = l2.*l3.*m3.*m4.*t4.*3.24e+2;
t186 = l3.*l4.*m2.*m4.*t3.*4.05e+2;
t187 = l3.*l4.*m3.*m4.*t3.*4.05e+2;
t197 = l1.*l2.*t3.*t15.*7.2e+1;
t198 = l3.*l4.*t3.*t17.*1.08e+2;
t199 = l3.*l4.*t3.*t15.*2.16e+2;
t200 = l2.*l3.*t4.*t14.*2.7e+2;
t201 = l2.*l3.*t4.*t15.*3.24e+2;
t202 = l1.*l4.*t3.*t17.*4.32e+2;
t224 = l2.*l3.*m2.*m4.*t4.*5.4e+2;
t225 = l1.*l4.*m1.*m3.*t3.*5.76e+2;
t226 = l1.*l4.*m1.*m4.*t3.*7.2e+2;
t227 = l2.*l3.*m2.*m3.*t4.*8.1e+2;
t228 = l1.*l4.*m2.*m3.*t3.*1.296e+3;
t229 = l1.*l4.*m2.*m4.*t3.*1.62e+3;
t230 = l1.*l4.*m3.*m4.*t3.*1.62e+3;
t283 = l1.*l4.*t3.*t15.*8.64e+2;
t40 = cos(t28);
t47 = cos(t41);
t48 = cos(t42);
t49 = cos(t43);
t50 = cos(t44);
t51 = cos(t45);
t52 = cos(t46);
t54 = t28+th2;
t58 = t39+th2;
t63 = cos(t57);
t65 = cos(t59);
t69 = t57+th4;
t73 = cos(t53);
t86 = t19+t39;
t94 = cos(t60);
t95 = cos(t61);
t96 = cos(t62);
t97 = t9+t28;
t98 = t21+t42;
t99 = t8+t28;
t105 = cos(t68);
t107 = cos(t70);
t135 = l2.*m2.*m3.*t26.*4.0;
t136 = l2.*m2.*m4.*t26.*5.0;
t137 = cos(t85);
t143 = t36+t45;
t144 = t38+t42;
t145 = -t100;
t146 = -t102;
t147 = -t104;
t157 = -t126;
t158 = -t127;
t159 = -t128;
t162 = l2.*m3.*m4.*t26.*4.5e+1;
t163 = -t139;
t164 = -t141;
t165 = t8+t53;
t169 = l2.*m2.*m4.*t24.*1.8e+1;
t170 = l2.*m3.*m4.*t24.*2.7e+1;
t171 = l2.*m3.*m4.*t23.*4.5e+1;
t189 = l2.*t17.*t26.*1.2e+1;
t190 = l2.*t15.*t26.*2.4e+1;
t191 = l1.*l3.*m1.*m4.*t24.*3.6e+1;
t192 = l1.*l3.*m3.*m4.*t24.*5.4e+1;
t193 = t22.*t71;
t194 = t22.*t72;
t195 = l2.*t15.*t23.*1.8e+1;
t196 = t23.*t75;
t206 = m1.*m2.*m4.*t24.*7.2e+1;
t207 = m1.*m3.*m4.*t24.*1.08e+2;
t208 = m1.*m2.*m4.*t24.*1.44e+2;
t209 = m1.*m3.*m4.*t23.*1.8e+2;
t210 = m1.*m3.*m4.*t24.*2.16e+2;
t211 = m2.*m3.*m4.*t23.*2.7e+2;
t212 = m2.*m3.*m4.*t24.*2.7e+2;
t213 = m1.*m3.*m4.*t23.*3.6e+2;
t214 = m2.*m3.*m4.*t22.*4.5e+2;
t215 = -t178;
t216 = -t179;
t217 = -t181;
t218 = -t182;
t219 = -t183;
t220 = -t184;
t221 = -t185;
t222 = -t186;
t223 = -t187;
t235 = l1.*l3.*m2.*m4.*t24.*1.08e+2;
t236 = l2.*l4.*m2.*m4.*t22.*2.7e+2;
t237 = l2.*l4.*m2.*m3.*t22.*4.32e+2;
t238 = t16.*t22.*-7.2e+1;
t239 = t16.*t22.*-1.44e+2;
t244 = m2.*m3.*t11.*t26.*1.2e+1;
t245 = m2.*m4.*t11.*t26.*1.5e+1;
t246 = m1.*m3.*t12.*t27.*3.6e+1;
t247 = m3.*m4.*t12.*t27.*5.4e+1;
t248 = l1.*l3.*m1.*m3.*t27.*1.44e+2;
t249 = l2.*l4.*m2.*m3.*t26.*1.44e+2;
t250 = l2.*l4.*m2.*m4.*t26.*1.8e+2;
t251 = l1.*l3.*m3.*m4.*t27.*2.16e+2;
t252 = l1.*l3.*m2.*m3.*t27.*3.24e+2;
t254 = l2.*t17.*t23.*-1.8e+1;
t255 = m1.*t15.*t23.*7.2e+1;
t256 = m3.*t14.*t22.*7.2e+1;
t257 = t23.*t108;
t258 = t22.*t109;
t259 = m4.*t15.*t24.*8.1e+1;
t260 = m2.*t15.*t23.*1.08e+2;
t261 = m2.*t17.*t22.*1.08e+2;
t262 = m2.*t17.*t23.*1.08e+2;
t263 = m4.*t14.*t22.*1.35e+2;
t264 = m4.*t14.*t24.*1.35e+2;
t265 = m1.*t15.*t23.*1.44e+2;
t266 = m3.*t14.*t22.*1.44e+2;
t267 = t23.*t112;
t268 = t22.*t113;
t269 = m4.*t15.*t24.*1.62e+2;
t270 = m2.*t15.*t22.*1.8e+2;
t271 = m2.*t15.*t23.*2.16e+2;
t272 = m2.*t17.*t22.*2.16e+2;
t273 = m2.*t17.*t23.*2.16e+2;
t274 = t22.*t116;
t275 = m4.*t14.*t22.*2.7e+2;
t276 = m4.*t14.*t24.*2.7e+2;
t277 = m2.*t15.*t22.*3.6e+2;
t278 = t22.*t122;
t279 = -t198;
t280 = -t199;
t281 = -t200;
t282 = -t201;
t293 = m2.*m3.*m4.*t23.*5.4e+2;
t294 = m2.*m3.*m4.*t24.*5.4e+2;
t295 = m2.*m3.*m4.*t22.*9.0e+2;
t296 = -t224;
t297 = -t227;
t308 = l2.*l4.*t14.*t22.*1.08e+2;
t309 = t22.*t126;
t310 = t22.*t128;
t312 = t22.*t133;
t317 = t12.*t15.*t27.*2.7e+1;
t318 = t11.*t17.*t26.*3.6e+1;
t319 = l1.*l3.*t15.*t27.*1.08e+2;
t320 = l2.*l4.*t17.*t26.*4.32e+2;
t329 = m2.*m3.*t12.*t27.*8.1e+1;
t330 = m3.*m4.*t11.*t26.*1.35e+2;
t331 = m1.*m4.*t12.*t27.*2.16e+2;
t332 = m2.*m4.*t12.*t27.*4.86e+2;
t338 = l1.*l3.*m1.*m4.*t27.*8.64e+2;
t339 = l2.*l4.*m3.*m4.*t26.*1.62e+3;
t340 = l1.*l3.*m2.*m4.*t27.*1.944e+3;
t343 = m1.*t17.*t23.*-7.2e+1;
t344 = m3.*t17.*t22.*-7.2e+1;
t353 = m1.*t17.*t23.*-1.44e+2;
t354 = m3.*t17.*t22.*-1.44e+2;
t360 = m4.*t15.*t22.*-2.25e+2;
t364 = m4.*t15.*t22.*-4.5e+2;
t410 = t11.*t15.*t26.*7.2e+1;
t412 = l2.*l4.*t15.*t26.*8.64e+2;
t64 = cos(t58);
t74 = cos(t54);
t106 = cos(t69);
t138 = cos(t86);
t150 = cos(t97);
t151 = cos(t98);
t152 = cos(t99);
t172 = cos(t143);
t173 = cos(t144);
t188 = cos(t165);
t203 = -t169;
t204 = -t170;
t205 = -t171;
t231 = l2.*m2.*m4.*t73.*3.0;
t232 = l2.*l3.*m2.*m3.*t40.*3.6e+1;
t233 = -t191;
t234 = -t192;
t240 = l3.*m2.*m4.*t49.*9.0;
t241 = l3.*m3.*m4.*t49.*9.0;
t242 = l3.*m3.*m4.*t47.*3.0e+1;
t243 = l1.*m3.*m4.*t47.*4.5e+1;
t253 = -t195;
t284 = -t206;
t285 = -t207;
t286 = -t208;
t287 = -t209;
t288 = -t210;
t289 = -t211;
t290 = -t212;
t291 = -t213;
t292 = -t214;
t298 = l1.*l2.*m3.*m4.*t51.*2.7e+1;
t299 = l2.*l3.*m2.*m4.*t49.*2.7e+1;
t300 = l2.*l3.*m3.*m4.*t49.*2.7e+1;
t301 = l1.*l2.*m2.*m3.*t48.*3.6e+1;
t302 = l1.*l2.*m1.*m4.*t51.*3.6e+1;
t303 = l1.*l2.*m2.*m4.*t48.*4.5e+1;
t305 = l2.*m3.*m4.*t73.*9.0;
t306 = l2.*l3.*m2.*m4.*t40.*2.16e+2;
t307 = l2.*l3.*m3.*m4.*t40.*2.16e+2;
t311 = -t235;
t313 = l3.*t15.*t47.*1.2e+1;
t314 = l3.*t17.*t47.*1.2e+1;
t315 = l1.*t15.*t47.*1.8e+1;
t316 = l1.*t17.*t47.*1.8e+1;
t321 = l2.*m2.*m3.*t63.*1.2e+1;
t324 = l2.*m2.*m4.*t63.*1.5e+1;
t326 = l2.*m3.*m4.*t63.*4.5e+1;
t327 = -t244;
t328 = -t245;
t333 = -t248;
t334 = -t249;
t335 = -t250;
t336 = -t251;
t337 = -t252;
t341 = -t255;
t342 = -t256;
t345 = -t259;
t346 = -t260;
t347 = -t261;
t348 = -t262;
t349 = -t263;
t350 = -t264;
t351 = -t265;
t352 = -t266;
t355 = -t269;
t356 = -t270;
t357 = -t271;
t358 = -t272;
t359 = -t273;
t361 = -t275;
t362 = -t276;
t363 = -t277;
t365 = -t293;
t366 = -t294;
t367 = -t295;
t368 = l1.*l2.*t17.*t48.*3.6e+1;
t369 = l2.*l3.*t15.*t47.*3.6e+1;
t370 = l2.*l3.*t17.*t47.*3.6e+1;
t371 = l2.*m3.*m4.*t96.*9.0;
t372 = m1.*m3.*m4.*t96.*3.6e+1;
t373 = m2.*m3.*m4.*t96.*5.4e+1;
t382 = l1.*l2.*m2.*m4.*t51.*8.1e+1;
t383 = l3.*l4.*m3.*m4.*t51.*8.1e+1;
t384 = l2.*l3.*m3.*m4.*t47.*9.0e+1;
t385 = l2.*l3.*m1.*m3.*t52.*1.08e+2;
t386 = l3.*l4.*m2.*m3.*t48.*1.08e+2;
t387 = l3.*l4.*m1.*m4.*t51.*1.08e+2;
t388 = l1.*l2.*m3.*m4.*t48.*1.35e+2;
t389 = l3.*l4.*m2.*m4.*t48.*1.35e+2;
t390 = l2.*l3.*m2.*m3.*t52.*1.62e+2;
t391 = l2.*l3.*m2.*m4.*t50.*1.62e+2;
t392 = l2.*l3.*m3.*m4.*t50.*1.62e+2;
t393 = l2.*l3.*m1.*m4.*t52.*2.16e+2;
t394 = l2.*l3.*m2.*m3.*t50.*2.43e+2;
t395 = l3.*l4.*m2.*m4.*t51.*2.43e+2;
t396 = l2.*l3.*m2.*m4.*t52.*3.24e+2;
t397 = l1.*l4.*m3.*m4.*t51.*3.24e+2;
t398 = l3.*l4.*m3.*m4.*t48.*4.05e+2;
t399 = l1.*l4.*m2.*m3.*t48.*4.32e+2;
t400 = l1.*l4.*m1.*m4.*t51.*4.32e+2;
t401 = l2.*l3.*t15.*t40.*1.08e+2;
t404 = l2.*t17.*t63.*1.2e+1;
t407 = l2.*t15.*t63.*2.4e+1;
t408 = -t317;
t409 = -t318;
t411 = -t320;
t415 = l1.*m2.*m4.*t49.*(2.7e+1./2.0);
t416 = l1.*m3.*m4.*t49.*(2.7e+1./2.0);
t418 = -t330;
t419 = -t338;
t420 = -t339;
t421 = -t340;
t422 = m2.*t15.*t94.*3.6e+1;
t423 = m2.*t17.*t94.*3.6e+1;
t427 = l1.*l2.*t15.*t48.*7.2e+1;
t428 = l2.*l3.*t14.*t50.*8.1e+1;
t429 = l3.*l4.*t17.*t48.*1.08e+2;
t430 = l2.*l3.*t15.*t50.*1.62e+2;
t431 = l3.*l4.*t15.*t48.*2.16e+2;
t432 = l1.*l4.*t17.*t48.*4.32e+2;
t433 = m2.*m3.*t11.*t63.*3.6e+1;
t434 = m2.*m4.*t11.*t63.*4.5e+1;
t435 = m1.*m3.*m4.*t96.*7.2e+1;
t436 = m2.*m3.*m4.*t95.*8.1e+1;
t437 = m2.*m3.*m4.*t94.*9.0e+1;
t438 = m2.*m3.*m4.*t96.*1.08e+2;
t439 = m2.*m3.*m4.*t95.*1.62e+2;
t440 = m2.*m3.*m4.*t94.*1.8e+2;
t447 = l1.*l3.*m1.*m3.*t65.*4.32e+2;
t449 = l2.*l4.*m2.*m3.*t63.*4.32e+2;
t451 = l1.*l4.*m2.*m4.*t48.*5.4e+2;
t452 = l1.*l4.*m2.*m4.*t51.*9.72e+2;
t453 = l1.*l4.*m3.*m4.*t48.*1.62e+3;
t457 = m2.*m4.*t11.*t73.*9.0;
t459 = m3.*m4.*t11.*t73.*2.7e+1;
t462 = l2.*l4.*m2.*m4.*t73.*1.08e+2;
t464 = l2.*l4.*m3.*m4.*t73.*3.24e+2;
t467 = -t410;
t468 = -t412;
t473 = t11.*t17.*t63.*3.6e+1;
t474 = m2.*t15.*t94.*7.2e+1;
t475 = m2.*t17.*t94.*7.2e+1;
t476 = m4.*t14.*t95.*8.1e+1;
t477 = m4.*t15.*t95.*8.1e+1;
t479 = l1.*l3.*t15.*t65.*3.24e+2;
t481 = l2.*l4.*t17.*t63.*4.32e+2;
t482 = l1.*l4.*t15.*t48.*8.64e+2;
t483 = m1.*m3.*t12.*t65.*1.08e+2;
t484 = m3.*m4.*t11.*t63.*1.35e+2;
t485 = m3.*m4.*t12.*t65.*1.62e+2;
t486 = m1.*m4.*t12.*t65.*2.16e+2;
t487 = m2.*m3.*t12.*t65.*2.43e+2;
t488 = m2.*m4.*t12.*t65.*4.86e+2;
t492 = l2.*l4.*m2.*m4.*t63.*5.4e+2;
t493 = l1.*l3.*m3.*m4.*t65.*6.48e+2;
t494 = l1.*l3.*m1.*m4.*t65.*8.64e+2;
t495 = l1.*l3.*m2.*m3.*t65.*9.72e+2;
t498 = l2.*l4.*m3.*m4.*t63.*1.62e+3;
t499 = l1.*l3.*m2.*m4.*t65.*1.944e+3;
t500 = l1.*l3.*t15.*t94.*3.6e+1;
t501 = t83.*t94;
t503 = l1.*l3.*m3.*m4.*t94.*9.0e+1;
t504 = l2.*l4.*m1.*m4.*t96.*1.08e+2;
t505 = l2.*l4.*m2.*m4.*t96.*1.62e+2;
t512 = l2.*l3.*m2.*m3.*t105.*1.08e+2;
t513 = l2.*l3.*m2.*m3.*t107.*1.08e+2;
t514 = l2.*l3.*m2.*m4.*t105.*2.16e+2;
t515 = l2.*l3.*m3.*m4.*t107.*2.16e+2;
t524 = t11.*t15.*t63.*7.2e+1;
t525 = t12.*t15.*t65.*8.1e+1;
t529 = l2.*l4.*t15.*t63.*8.64e+2;
t536 = l1.*l3.*t17.*t94.*-3.6e+1;
t538 = l2.*l3.*m2.*m4.*t137.*1.62e+2;
t539 = l2.*l3.*m3.*m4.*t137.*1.62e+2;
t540 = l2.*l3.*m2.*m3.*t137.*2.43e+2;
t542 = l2.*l3.*t15.*t107.*1.08e+2;
t543 = l2.*l3.*t15.*t105.*3.24e+2;
t551 = l2.*l3.*m2.*m4.*t107.*6.48e+2;
t552 = l2.*l3.*m3.*m4.*t105.*6.48e+2;
t563 = m4.*t14.*t95.*(8.1e+1./2.0);
t564 = m4.*t15.*t95.*(8.1e+1./2.0);
t567 = l2.*l3.*t14.*t137.*8.1e+1;
t568 = l2.*l3.*t15.*t137.*1.62e+2;
t304 = -t231;
t322 = l3.*m2.*m4.*t64.*9.0;
t323 = l3.*m3.*m4.*t64.*9.0;
t325 = -t243;
t374 = -t298;
t375 = l2.*l3.*m2.*m4.*t64.*2.7e+1;
t376 = -t299;
t377 = l2.*l3.*m3.*m4.*t64.*2.7e+1;
t378 = -t300;
t379 = -t301;
t380 = -t302;
t381 = -t303;
t402 = -t305;
t403 = l3.*m3.*m4.*t150.*6.0;
t405 = -t315;
t406 = -t316;
t413 = -t321;
t414 = -t324;
t417 = -t326;
t424 = -t368;
t425 = -t369;
t426 = -t370;
t443 = -t382;
t444 = -t384;
t445 = -t388;
t446 = -t397;
t448 = -t399;
t450 = -t400;
t454 = -t401;
t456 = l1.*m3.*m4.*t150.*9.0;
t458 = m2.*m3.*t12.*t74.*2.7e+1;
t460 = m3.*m4.*t12.*t74.*5.4e+1;
t461 = l1.*l3.*m2.*m3.*t74.*1.08e+2;
t463 = l1.*l3.*m3.*m4.*t74.*2.16e+2;
t465 = -t404;
t466 = -t407;
t469 = l1.*m2.*m4.*t64.*(2.7e+1./2.0);
t470 = -t415;
t471 = l1.*m3.*m4.*t64.*(2.7e+1./2.0);
t472 = -t416;
t478 = -t427;
t480 = -t432;
t489 = m2.*m3.*m4.*t138.*8.1e+1;
t490 = m2.*m3.*m4.*t138.*1.62e+2;
t491 = -t451;
t496 = -t452;
t497 = -t453;
t502 = t12.*t15.*t74.*2.7e+1;
t506 = l1.*l3.*t15.*t74.*1.08e+2;
t507 = l2.*m2.*m4.*t172.*9.0;
t508 = l2.*m3.*m4.*t172.*9.0;
t511 = m2.*m4.*t12.*t74.*1.62e+2;
t516 = l2.*l3.*m2.*m3.*t106.*3.24e+2;
t517 = l1.*l3.*m2.*m4.*t74.*6.48e+2;
t520 = l2.*l3.*m3.*m4.*t150.*1.8e+1;
t521 = l1.*l2.*m2.*m4.*t151.*2.7e+1;
t522 = l1.*l2.*m3.*m4.*t151.*2.7e+1;
t523 = l2.*l3.*m2.*m3.*t152.*5.4e+1;
t526 = m4.*t14.*t138.*8.1e+1;
t527 = m4.*t15.*t138.*8.1e+1;
t528 = -t482;
t530 = -t483;
t531 = -t485;
t532 = -t486;
t533 = -t487;
t534 = -t488;
t535 = -t500;
t537 = -t503;
t544 = l2.*l3.*t15.*t106.*3.24e+2;
t546 = -t512;
t547 = -t513;
t548 = -t514;
t549 = -t515;
t550 = l2.*l3.*m2.*m4.*t106.*6.48e+2;
t553 = l2.*l3.*m3.*m4.*t106.*6.48e+2;
t554 = m2.*m3.*m4.*t188.*1.8e+1;
t555 = m2.*m3.*m4.*t188.*3.6e+1;
t557 = l3.*l4.*m2.*m4.*t151.*8.1e+1;
t558 = l3.*l4.*m3.*m4.*t151.*8.1e+1;
t559 = l2.*l3.*m2.*m4.*t152.*1.08e+2;
t560 = l1.*l4.*m2.*m4.*t151.*3.24e+2;
t561 = l1.*l4.*m3.*m4.*t151.*3.24e+2;
t562 = -t525;
t565 = l1.*l3.*m3.*m4.*t188.*1.8e+1;
t566 = l2.*l4.*m2.*m4.*t188.*5.4e+1;
t569 = -t543;
t570 = -t551;
t571 = -t552;
t572 = m2.*m4.*t11.*t172.*2.7e+1;
t573 = m3.*m4.*t11.*t172.*2.7e+1;
t579 = l1.*l3.*m2.*m3.*t173.*3.24e+2;
t580 = l2.*l4.*m2.*m4.*t172.*3.24e+2;
t581 = l2.*l4.*m3.*m4.*t172.*3.24e+2;
t582 = m4.*t14.*t138.*(8.1e+1./2.0);
t583 = m4.*t15.*t138.*(8.1e+1./2.0);
t585 = l1.*l3.*t15.*t173.*3.24e+2;
t588 = m2.*m3.*t12.*t173.*8.1e+1;
t589 = m2.*m4.*t12.*t173.*1.62e+2;
t590 = m3.*m4.*t12.*t173.*1.62e+2;
t594 = l1.*l3.*m2.*m4.*t173.*6.48e+2;
t595 = l1.*l3.*m3.*m4.*t173.*6.48e+2;
t596 = t12.*t15.*t173.*8.1e+1;
t441 = -t375;
t442 = -t377;
t455 = -t403;
t509 = -t458;
t510 = -t460;
t518 = -t469;
t519 = -t471;
t541 = -t506;
t545 = -t511;
t556 = -t523;
t574 = -t554;
t575 = -t555;
t576 = -t557;
t577 = -t558;
t578 = -t559;
t584 = -t566;
t586 = -t572;
t587 = -t573;
t591 = -t579;
t592 = -t580;
t593 = -t581;
t597 = -t585;
t598 = -t594;
t599 = -t595;
t600 = t55+t56+t75+t76+t77+t101+t103+t140+t142+t148+t203+t204+t205+t253+t254+t325+t371+t405+t406+t456+t470+t472+t518+t519;
t601 = t135+t136+t145+t146+t147+t162+t163+t164+t189+t190+t240+t241+t242+t304+t313+t314+t322+t323+t402+t413+t414+t417+t455+t465+t466+t507+t508;
t602 = t71+t78+t79+t82+t108+t109+t110+t111+t114+t115+t116+t119+t124+t238+t284+t285+t287+t289+t290+t292+t341+t342+t343+t344+t345+t346+t347+t348+t349+t350+t356+t360+t372+t373+t422+t423+t436+t437+t489+t563+t564+t574+t582+t583;
t603 = t72+t80+t81+t112+t113+t117+t118+t120+t121+t122+t123+t125+t149+t239+t286+t288+t291+t351+t352+t353+t354+t355+t357+t358+t359+t361+t362+t363+t364+t365+t366+t367+t435+t438+t439+t440+t474+t475+t476+t477+t490+t526+t527+t575;
t606 = t66+t67+t83+t84+t87+t89+t90+t153+t154+t155+t166+t167+t168+t174+t175+t176+t177+t180+t197+t233+t234+t311+t327+t328+t374+t376+t378+t379+t380+t381+t409+t418+t424+t425+t426+t433+t434+t441+t442+t443+t444+t445+t457+t459+t467+t473+t478+t484+t520+t521+t522+t524+t535+t536+t537+t565+t586+t587;
t611 = t202+t225+t226+t228+t229+t230+t232+t283+t306+t307+t319+t333+t334+t335+t336+t337+t411+t419+t420+t421+t446+t447+t448+t449+t450+t454+t461+t462+t463+t464+t468+t479+t480+t481+t491+t492+t493+t494+t495+t496+t497+t498+t499+t516+t517+t528+t529+t541+t542+t544+t546+t547+t548+t549+t550+t553+t560+t561+t569+t570+t571+t591+t592+t593+t597+t598+t599;
t612 = t129+t130+t131+t132+t157+t158+t159+t160+t161+t215+t216+t217+t218+t219+t220+t221+t222+t223+t236+t237+t246+t247+t279+t280+t281+t282+t296+t297+t308+t309+t310+t312+t329+t331+t332+t383+t385+t386+t387+t389+t390+t391+t392+t393+t394+t395+t396+t398+t408+t428+t429+t430+t431+t502+t504+t505+t509+t510+t530+t531+t532+t533+t534+t538+t539+t540+t545+t556+t562+t567+t568+t576+t577+t578+t584+t588+t589+t590+t596;
t604 = 1.0./t602;
t605 = 1.0./t603;
t607 = t30.*t31.*t600.*t604.*1.2e+1;
t609 = t29.*t31.*t33.*t601.*t605.*3.6e+1;
t613 = t29.*t32.*t33.*t605.*t606.*1.2e+1;
t615 = t29.*t31.*t33.*t35.*t605.*t611;
t616 = t31.*t34.*t35.*t605.*t612.*4.0;
t608 = -t607;
t610 = -t609;
t614 = -t613;

out = zeros(size(x,1), size(u,1), size(x,2));

out(5,1,:) = t30.*t605.*(t15.*3.0e+1+t17.*1.8e+1+m2.*m3.*1.6e+1+m2.*m4.*3.0e+1+m3.*m4.*7.5e+1-t15.*t23.*1.8e+1-t17.*t23.*1.8e+1-m2.*m4.*t24.*1.8e+1-m3.*m4.*t23.*4.5e+1-m3.*m4.*t24.*2.7e+1+m3.*m4.*t96.*9.0).*2.4e+1;
out(6,1,:) = t608;
out(7,1,:) = t610;
out(8,1,:) = -t29.*t33.*t35.*t605.*(l4.*t5.*t6.*t15.*1.728e+3+l4.*t5.*t6.*t17.*8.64e+2+l4.*m2.*m3.*t2.*t3.*2.88e+2+l4.*m2.*m4.*t2.*t3.*3.6e+2+l4.*m2.*m3.*t5.*t6.*5.76e+2+l4.*m2.*m4.*t5.*t6.*7.2e+2+l4.*m3.*m4.*t5.*t6.*3.24e+3-l3.*t3.*t5.*t7.*t15.*4.32e+2+l3.*t4.*t5.*t6.*t15.*8.64e+2+l3.*m2.*m3.*t2.*t3.*t4.*1.44e+2+l3.*m2.*m3.*t2.*t6.*t7.*2.88e+2-l3.*m2.*m3.*t3.*t5.*t7.*5.76e+2+l3.*m2.*m3.*t4.*t5.*t6.*2.88e+2+l3.*m2.*m4.*t2.*t6.*t7.*8.64e+2-l3.*m2.*m4.*t3.*t5.*t7.*1.728e+3-l3.*m3.*m4.*t3.*t5.*t7.*1.728e+3+l3.*m3.*m4.*t4.*t5.*t6.*8.64e+2-l4.*m2.*m4.*t2.*t3.*t24.*2.16e+2+l4.*m2.*m4.*t2.*t6.*t25.*2.16e+2-l4.*m2.*m4.*t3.*t5.*t25.*4.32e+2-l4.*m3.*m4.*t3.*t5.*t25.*6.48e+2-l4.*m2.*m4.*t5.*t6.*t24.*4.32e+2-l4.*m3.*m4.*t5.*t6.*t24.*6.48e+2);

out(5,2,:) = t608;
out(6,2,:) = t30.*t32.*t605.*(t156+t10.*t15.*3.0e+1+t11.*t15.*3.0e+1+t10.*t17.*1.8e+1+m1.*m3.*t10.*1.6e+1+m1.*m4.*t10.*3.0e+1+m2.*m3.*t10.*4.8e+1+m2.*m3.*t11.*1.6e+1+m2.*m4.*t10.*9.0e+1+m2.*m4.*t11.*3.0e+1+m3.*m4.*t10.*7.5e+1+m3.*m4.*t11.*7.5e+1-t11.*t15.*t23.*1.8e+1-t11.*t17.*t23.*1.8e+1-t10.*t15.*t94.*1.8e+1-t10.*t17.*t94.*1.8e+1+l1.*l2.*t2.*t15.*6.0e+1+l1.*l2.*t2.*t17.*3.6e+1-l1.*l2.*t15.*t47.*3.6e+1-l1.*l2.*t17.*t47.*3.6e+1-m1.*m4.*t10.*t24.*1.8e+1-m2.*m4.*t10.*t24.*5.4e+1-m2.*m4.*t11.*t24.*1.8e+1-m3.*m4.*t10.*t24.*2.7e+1-m3.*m4.*t11.*t23.*4.5e+1-m3.*m4.*t11.*t24.*2.7e+1-m3.*m4.*t10.*t94.*4.5e+1+m3.*m4.*t11.*t96.*9.0+m3.*m4.*t10.*t188.*9.0+l1.*l2.*m2.*m3.*t2.*4.8e+1+l1.*l2.*m2.*m4.*t2.*9.0e+1+l1.*l2.*m3.*m4.*t2.*1.5e+2-l1.*l2.*m3.*m4.*t47.*9.0e+1-l1.*l2.*m2.*m4.*t49.*2.7e+1-l1.*l2.*m3.*m4.*t49.*2.7e+1-l1.*l2.*m2.*m4.*t64.*2.7e+1-l1.*l2.*m3.*m4.*t64.*2.7e+1+l1.*l2.*m3.*m4.*t150.*1.8e+1).*2.4e+1;
out(7,2,:) = t614;
out(8,2,:) = t615;

out(5,3,:) = t610;
out(6,3,:) = t614;
out(7,3,:) = t32.*t34.*t605.*(t156+t11.*t14.*3.0e+1+t11.*t15.*7.2e+1+t12.*t15.*3.0e+1+t12.*t17.*1.8e+1+m1.*m2.*t11.*1.6e+1+m1.*m3.*t11.*4.8e+1+m1.*m3.*t12.*1.6e+1+m1.*m4.*t11.*3.0e+1+m2.*m3.*t11.*1.2e+2+m1.*m4.*t12.*3.0e+1+m2.*m3.*t12.*4.8e+1+m2.*m4.*t11.*7.5e+1+m2.*m4.*t12.*9.0e+1+m3.*m4.*t11.*9.0e+1+m3.*m4.*t12.*7.5e+1-t11.*t14.*t22.*1.8e+1-t11.*t15.*t22.*7.2e+1-t11.*t17.*t22.*1.8e+1-t12.*t15.*t94.*1.8e+1-t12.*t17.*t94.*1.8e+1+l2.*l3.*t3.*t15.*7.2e+1+l2.*l3.*t3.*t17.*3.6e+1-l2.*l3.*t15.*t48.*7.2e+1-l2.*l3.*t17.*t48.*3.6e+1-m2.*m3.*t11.*t22.*7.2e+1-m2.*m4.*t11.*t22.*4.5e+1-m3.*m4.*t11.*t22.*9.0e+1-m1.*m4.*t12.*t24.*1.8e+1-m2.*m4.*t12.*t24.*5.4e+1-m3.*m4.*t12.*t24.*2.7e+1-m1.*m4.*t11.*t96.*1.8e+1-m2.*m4.*t11.*t96.*2.7e+1-m3.*m4.*t12.*t94.*4.5e+1+m2.*m4.*t11.*t188.*9.0+m3.*m4.*t12.*t188.*9.0+l2.*l3.*m1.*m3.*t3.*4.8e+1+l2.*l3.*m1.*m4.*t3.*6.0e+1+l2.*l3.*m2.*m3.*t3.*1.08e+2+l2.*l3.*m2.*m4.*t3.*1.35e+2+l2.*l3.*m3.*m4.*t3.*1.35e+2-l2.*l3.*m2.*m3.*t48.*3.6e+1-l2.*l3.*m2.*m4.*t48.*4.5e+1-l2.*l3.*m3.*m4.*t48.*1.35e+2-l2.*l3.*m1.*m4.*t51.*3.6e+1-l2.*l3.*m2.*m4.*t51.*8.1e+1-l2.*l3.*m3.*m4.*t51.*2.7e+1+l2.*l3.*m2.*m4.*t151.*2.7e+1+l2.*l3.*m3.*m4.*t151.*2.7e+1).*2.4e+1;
out(8,3,:) = t616;

out(5,4,:) = t29.*t33.*t35.*t605.*(l4.*t5.*t6.*t15.*2.4e+1+l4.*t5.*t6.*t17.*1.2e+1+l4.*m2.*m3.*t2.*t3.*4.0+l4.*m2.*m4.*t2.*t3.*5.0+l4.*m2.*m3.*t5.*t6.*8.0+l4.*m2.*m4.*t5.*t6.*1.0e+1+l4.*m3.*m4.*t5.*t6.*4.5e+1-l3.*t3.*t5.*t7.*t15.*6.0+l3.*t4.*t5.*t6.*t15.*1.2e+1+l3.*m2.*m3.*t2.*t3.*t4.*2.0+l3.*m2.*m3.*t2.*t6.*t7.*4.0-l3.*m2.*m3.*t3.*t5.*t7.*8.0+l3.*m2.*m3.*t4.*t5.*t6.*4.0+l3.*m2.*m4.*t2.*t6.*t7.*1.2e+1-l3.*m2.*m4.*t3.*t5.*t7.*2.4e+1-l3.*m3.*m4.*t3.*t5.*t7.*2.4e+1+l3.*m3.*m4.*t4.*t5.*t6.*1.2e+1-l4.*m2.*m4.*t2.*t3.*t24.*3.0+l4.*m2.*m4.*t2.*t6.*t25.*3.0-l4.*m2.*m4.*t3.*t5.*t25.*6.0-l4.*m3.*m4.*t3.*t5.*t25.*9.0-l4.*m2.*m4.*t5.*t6.*t24.*6.0-l4.*m3.*m4.*t5.*t6.*t24.*9.0).*-7.2e+1;
out(6,4,:) = t615;
out(7,4,:) = t616;
out(8,4,:) = (t34.*t605.*(t12.*t16.*1.8e+1+t13.*t18.*1.8e+1+t12.*t108+t12.*t109+t12.*t114+t12.*t343+t12.*t344+t12.*t423+m1.*t12.*t15.*3.0e+1+m2.*t12.*t15.*7.5e+1+m3.*t12.*t14.*3.0e+1+m4.*t12.*t14.*9.0e+1+m1.*t13.*t17.*3.0e+1+m4.*t12.*t15.*9.0e+1+m4.*t13.*t14.*3.0e+1+m2.*t13.*t17.*7.5e+1+m4.*t13.*t15.*7.2e+1+m3.*t13.*t17.*9.0e+1-t12.*t16.*t22.*1.8e+1-t13.*t18.*t22.*1.8e+1+m1.*m2.*m3.*t12.*1.6e+1+m1.*m2.*m4.*t12.*4.8e+1+m1.*m2.*m4.*t13.*1.6e+1+m1.*m3.*m4.*t12.*1.2e+2+m1.*m3.*m4.*t13.*4.8e+1+m2.*m3.*m4.*t12.*3.0e+2+m2.*m3.*m4.*t13.*1.2e+2+l3.*l4.*t4.*t108+l3.*l4.*t4.*t114-m1.*t12.*t15.*t23.*1.8e+1-m2.*t12.*t15.*t22.*4.5e+1-m3.*t12.*t14.*t22.*1.8e+1-m2.*t12.*t15.*t23.*2.7e+1-m4.*t12.*t14.*t22.*5.4e+1-m2.*t12.*t17.*t22.*1.08e+2-m4.*t12.*t15.*t22.*9.0e+1-m4.*t13.*t14.*t22.*1.8e+1-m2.*t12.*t17.*t23.*1.08e+2-m2.*t13.*t17.*t22.*4.5e+1-m4.*t13.*t15.*t22.*7.2e+1-m3.*t13.*t17.*t22.*9.0e+1+m2.*t12.*t15.*t94.*9.0-m1.*t13.*t17.*t96.*1.8e+1-m2.*t13.*t17.*t96.*2.7e+1+m2.*t13.*t17.*t188.*9.0+l3.*l4.*m4.*t4.*t14.*9.0e+1+l3.*l4.*m4.*t4.*t15.*1.08e+2+l3.*l4.*m3.*t4.*t17.*1.08e+2-l3.*l4.*m4.*t14.*t50.*2.7e+1-l3.*l4.*m2.*t17.*t50.*5.4e+1-l3.*l4.*m4.*t15.*t50.*5.4e+1-l3.*l4.*m1.*t17.*t52.*7.2e+1-l3.*l4.*m3.*t17.*t50.*5.4e+1-l3.*l4.*m2.*t17.*t52.*1.08e+2-l3.*l4.*m4.*t14.*t137.*2.7e+1-l3.*l4.*m2.*t17.*t137.*5.4e+1-l3.*l4.*m4.*t15.*t137.*5.4e+1-l3.*l4.*m3.*t17.*t137.*5.4e+1+l3.*l4.*m2.*t17.*t152.*3.6e+1-m1.*m3.*m4.*t12.*t23.*7.2e+1-m2.*m3.*m4.*t12.*t22.*1.8e+2-m2.*m3.*m4.*t12.*t23.*1.08e+2-m2.*m3.*m4.*t13.*t22.*7.2e+1+m2.*m3.*m4.*t12.*t94.*3.6e+1+l3.*l4.*m1.*m2.*m4.*t4.*4.8e+1+l3.*l4.*m1.*m3.*m4.*t4.*1.08e+2+l3.*l4.*m2.*m3.*m4.*t4.*2.7e+2-l3.*l4.*m2.*m3.*m4.*t50.*8.1e+1-l3.*l4.*m1.*m3.*m4.*t52.*3.6e+1-l3.*l4.*m2.*m3.*m4.*t52.*5.4e+1-l3.*l4.*m2.*m3.*m4.*t137.*8.1e+1+l3.*l4.*m2.*m3.*m4.*t152.*1.8e+1).*2.4e+1)./(m4.*t13);

end
