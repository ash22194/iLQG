/*
 * Example of how to use the mxGPUArray API in a MEX file.  This example shows
 * how to write a MEX function that takes a gpuArray input and returns a
 * gpuArray output, e.g. B=mexFunction(A).
 *
 * Copyright 2012 The MathWorks, Inc.
 */

#include "mex.h"
#include "matrix.h"
#include "gpu/mxGPUArray.h"
#include <iostream>

 /*
  * Device code
  */

__global__ void dyn0_mex_continuous(double* const d_PH, // outputs
    double* const PH, double* const TH, double* const OM, double* const dPH, double* const dTH,
    double* const dPS, double* const dOM, double* const dNE, // input states
    double const* const in1, double const* const in2, // input actions
    const double mw, const double Iwxx, const double Iwyy, const double Iwzz, 
    const double mf, const double Ifxx, const double Ifyy, const double Ifzz, 
    const double mt, const double Itxx, const double Ityy, const double Itzz, 
    const double nw, const double nt, const double rw, const double rf, const double rt, const double fcoeff, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3]
                        * grid_size[4] * grid_size[5] * grid_size[6] * grid_size[7];
    int index4 = (grid_size[3] > 1) ? index : 0;

    while (index < num_elements)
    {
        d_PH[index] = dPH[index4];
        
        index = index + num_threads;
        index4 = (grid_size[3] > 1) ? index : 0;
    }
}

__global__ void dyn1_mex_continuous(double* const d_TH, // outputs
    double* const PH, double* const TH, double* const OM, double* const dPH, double* const dTH,
    double* const dPS, double* const dOM, double* const dNE, // input states
    double const* const in1, double const* const in2, // input actions
    const double mw, const double Iwxx, const double Iwyy, const double Iwzz, 
    const double mf, const double Ifxx, const double Ifyy, const double Ifzz, 
    const double mt, const double Itxx, const double Ityy, const double Itzz, 
    const double nw, const double nt, const double rw, const double rf, const double rt, const double fcoeff, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3]
                        * grid_size[4] * grid_size[5] * grid_size[6] * grid_size[7];
    int index5 = (grid_size[4] > 1) ? index : 0;

    while (index < num_elements)
    {
        d_TH[index] = dTH[index5];
        
        index = index + num_threads;
        index5 = (grid_size[4] > 1) ? index : 0;
    }
}

__global__ void dyn2_mex_continuous(double* const d_OM, // outputs
    double* const PH, double* const TH, double* const OM, double* const dPH, double* const dTH,
    double* const dPS, double* const dOM, double* const dNE, // input states
    double const* const in1, double const* const in2, // input actions
    const double mw, const double Iwxx, const double Iwyy, const double Iwzz, 
    const double mf, const double Ifxx, const double Ifyy, const double Ifzz, 
    const double mt, const double Itxx, const double Ityy, const double Itzz, 
    const double nw, const double nt, const double rw, const double rf, const double rt, const double fcoeff, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3]
                        * grid_size[4] * grid_size[5] * grid_size[6] * grid_size[7];
    int index7 = (grid_size[6] > 1) ? index : 0;

    while (index < num_elements)
    {
        d_OM[index] = dOM[index7];
        
        index = index + num_threads;
        index7 = (grid_size[6] > 1) ? index : 0;
    }
}

__global__ void dyn3_mex_continuous(double* const d_dPH, // outputs
    double* const PH, double* const TH, double* const OM, double* const dPH, double* const dTH,
    double* const dPS, double* const dOM, double* const dNE, // input states
    double const* const in1, double const* const in2, // input actions
    const double mw, const double Iwxx, const double Iwyy, const double Iwzz, 
    const double mf, const double Ifxx, const double Ifyy, const double Ifzz, 
    const double mt, const double Itxx, const double Ityy, const double Itzz, 
    const double nw, const double nt, const double rw, const double rf, const double rt, const double fcoeff, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3]
                        * grid_size[4] * grid_size[5] * grid_size[6] * grid_size[7];
    double th, ps, om, ne, dph, dth, dps, dom, dne, Tneta, Tomega;
    ps = 0; ne = 0;

    // int index1 = (grid_size[0] > 1) ? index : 0;
    int index2 = (grid_size[1] > 1) ? index : 0;
    int index3 = (grid_size[2] > 1) ? index : 0;
    int index4 = (grid_size[3] > 1) ? index : 0;
    int index5 = (grid_size[4] > 1) ? index : 0;
    int index6 = (grid_size[5] > 1) ? index : 0;
    int index7 = (grid_size[6] > 1) ? index : 0;
    int index8 = (grid_size[7] > 1) ? index : 0;

    int uindex1 = (active_actions[0] == 1) ? index : 0;
    int uindex2 = (active_actions[1] == 1) ? index : 0;

    while (index < num_elements)
    {
        th = TH[index2];
        om = OM[index3];
        dph = dPH[index4];
        dth = dTH[index5];
        dps = dPS[index6];
        dom = dOM[index7];
        dne = dNE[index8];

        Tomega = nw*in1[uindex1];
        Tneta = nt*in2[uindex2];

        double t2 = cos(ph);
        double t3 = cos(th);
        double t4 = cos(ps);
        double t5 = cos(om);
        double t6 = cos(ne);
        double t7 = sin(ph);
        double t8 = sin(th);
        double t9 = sin(ps);
        double t10 = sin(om);
        double t11 = sin(ne);
        double t12 = mf+mt;
        double t13 = mf*rf;
        double t14 = mt*rt;
        double t15 = rf+rt;
        double t16 = Ifxx*2.0;
        double t17 = Ifxx*4.0;
        double t18 = Ifyy*2.0;
        double t19 = Ifyy*4.0;
        double t20 = Ifzz*2.0;
        double t21 = Ifzz*4.0;
        double t22 = Itxx*2.0;
        double t23 = Itxx*3.0;
        double t24 = Itxx*4.0;
        double t25 = Ityy*2.0;
        double t26 = Ityy*3.0;
        double t27 = Ityy*4.0;
        double t28 = Itzz*2.0;
        double t29 = Itzz*4.0;
        double t30 = Itzz*Itzz;
        double t31 = Iwxx*2.0;
        double t32 = Iwxx*4.0;
        double t33 = Iwyy*2.0;
        double t34 = Iwyy*4.0;
        double t35 = Iwzz*2.0;
        double t36 = Iwzz*4.0;
        double t37 = fcoeff*2.0;

        double t40 = mt*2.0;
        double t42 = rf*rf;

        double t44 = rw*rw;
        double t45 = Tomega*4.0;
        double t46 = th*2.0;
        double t47 = ps*2.0;
        double t48 = om*2.0;
        double t49 = ne*2.0;
        double t50 = dom*2.0;
        double t51 = dph*dph;
        double t52 = dth*dth;
        double t53 = dom*dom;
        double t79 = Ifyy*8.0;
        double t80 = -Ifzz;
        double t83 = -Ityy;
        double t87 = -Itzz;
        double t92 = Iwyy*8.0;
        double t93 = -Iwzz;
        double t109 = Ifxx/2.0;
        double t110 = Ifzz/2.0;
        double t111 = Itxx/4.0;
        double t112 = Ityy/4.0;
        double t113 = Itzz/2.0;
        double t114 = Iwxx/2.0;
        double t115 = Iwzz/2.0;
        double t54 = cos(t46);
        double t55 = cos(t47);
        double t56 = cos(t48);
        double t57 = cos(t49);
        double t58 = t2*t2;
        double t59 = t3*t3;
        double t60 = t4*t4;
        double t61 = t5*t5;
        double t62 = t5*t5*t5;
        double t63 = t6*t6;
        double t64 = t13*2.0;
        double t65 = t13*4.0;
        double t66 = t14*2.0;
        double t67 = sin(t46);
        double t68 = sin(t47);
        double t69 = sin(t48);
        double t70 = sin(t49);
        double t71 = t7*t7;
        double t72 = t8*t8;
        double t73 = t9*t9;
        double t74 = t10*t10;
        double t75 = t11*t11;
        double t76 = mw+t12;
        double t77 = -t16;
        double t78 = -t17;
        double t81 = -t20;
        double t82 = -t21;
        double t84 = -t25;
        double t85 = -t26;
        double t86 = -t27;
        double t88 = -t28;
        double t89 = -t29;
        double t90 = -t31;
        double t91 = -t32;
        double t94 = -t35;
        double t95 = -t36;
        double t96 = t8*dph;
        double t97 = -t45;
        double t98 = rf*t12;
        double t99 = mt*t15;
        double t100 = t2*t5;
        double t101 = t2*t10;
        double t102 = t5*t7;
        double t103 = t6*t8;
        double t104 = t2*dph*dps;
        double t105 = t7*t10;
        double t106 = t8*t11;
        double t107 = rf*t13;
        double t108 = t15*t15;
        double t116 = rf*t14*4.0;
        double t117 = rf*t14*6.0;
        double t118 = Ifyy*Iwyy*t8;
        double t120 = t15*t40;
        double t122 = -t52;
        double t124 = Itxx+t83;
        double t125 = Iwxx+t93;
        double t133 = rt*t14*3.0;
        double t134 = rt*t14*4.0;
        double t135 = rf*t14*8.0;
        double t140 = t3*t6*t10;
        double t141 = t51+t52;
        double t144 = t3*t10*t11;
        double t148 = t20+t28;
        double t155 = t12*2.0;
        double t156 = t12*3.0;
        double t157 = t2*t3*dph*dth*2.0;
        double t158 = t7*t8*t52;
        double t159 = t3*t7*dph*dth*2.0;
        double t179 = Iwyy*mt*t8*t42;
        double t180 = Iwyy*rt*t8*t14;
        double t181 = (Itxx*Iwyy*t8)/2.0;
        double t182 = (Ityy*Iwyy*t8)/2.0;
        double t183 = rf*t8*t14*t33;
        double t191 = Iwyy*rw*t3*t7*t8;
        double t219 = t12*t42*4.0;
        double t119 = t98*2.0;
        double t121 = t99*4.0;
        double t123 = t50*t96;
        double t126 = rw*t76;
        double t127 = t96+dps;
        double t128 = t96+dom;
        double t129 = t54*3.0;
        double t130 = rf*t64;
        double t131 = rf*t65;
        double t132 = rt*t66;

        double t137 = Itxx*t63;
        double t138 = Ityy*t75;
        double t139 = t8*t100;
        double t142 = t8*t101;
        double t143 = t8*t102;
        double t145 = t8*t105;
        double t146 = t15*t99;
        double t149 = t22*t63;
        double t150 = t27*t63;
        double t151 = t35*t60;
        double t152 = t24*t75;
        double t153 = t25*t75;
        double t154 = t31*t73;
        double t160 = Iwyy*t13*t101;
        double t161 = Iwyy*mt*rf*t101;
        double t162 = Iwyy*t14*t101;
        double t163 = Iwyy*t13*t105;
        double t164 = Iwyy*mt*rf*t105;
        double t165 = Iwyy*t14*t105;
        double t166 = t40*t108;
        double t168 = t51*t54;
        double t169 = t51*t59;
        double t171 = t107/2.0;
        double t172 = t22+t84;
        double t173 = t23+t85;
        double t174 = t24+t86;
        double t175 = t124*t124;
        double t176 = t31+t94;
        double t177 = t32+t95;
        double t178 = Iwyy*t8*t107;
        double t184 = t50+t96;
        double t186 = -t140;
        double t187 = -t157;
        double t189 = t14+t98;
        double t190 = t13+t99;
        double t203 = t56*t59*2.0;
        double t205 = t42*t155;
        double t206 = t42*t156;

        double t208 = Ityy*Iwyy*t11*t140;
        double t210 = t64+t120;
        double t212 = t57*t124;
        double t213 = t55*t125;
        double t214 = t68*t125;
        double t215 = t44*t58*t76;

        double t217 = t57*t182;
        double t218 = t44*t71*t76;
        double t223 = t2*t8*t141;
        double t225 = t61*t148;
        double t226 = t10*t57*t67*4.0;
        double t232 = Itxx*Iwyy*t8*t57*(-1.0/2.0);
        double t236 = t5*t6*t11*t124;
        double t244 = t103+t144;
        double t246 = t8*t70*t124;
        double t253 = t6*t11*t53*t124;
        double t254 = t4*t9*t52*t125;

        double t275 = t3*t10*t70*t124;
        double t287 = t5*t67*t70*t124;
        double t288 = t3*t69*t70*t124;

        double t333 = t6*t11*t61*t122*t124;
        double t336 = t3*t6*t11*t61*t87*t124;
        double t147 = -t129;
        double t167 = t15*t121;

        double t185 = -t139;
        double t188 = -t145;

        double t195 = Iwyy*t13*t143;
        double t196 = Iwyy*mt*rf*t143;
        double t197 = Iwyy*t14*t143;
        double t198 = t2*t8*t126;
        double t199 = t7*t8*t126;
        double t200 = t7*t127*dph;
        double t201 = t146/2.0;
        double t202 = -t168;
        double t204 = t190*t190;
        double t209 = t66+t119;
        double t211 = t65+t121;
        double t227 = t5*t189;
        double t228 = t5*t190;
        double t229 = t212*2.0;
        double t230 = t213*2.0;
        double t231 = t213*4.0;
        double t233 = t51+t53+t123;
        double t234 = t57*t173;
        double t235 = t55*t176;
        double t237 = t68*t176;
        double t239 = t107+t146;
        double t240 = Itxx*Iwyy*t11*t186;
        double t241 = Itzz+t212;
        double t242 = t101+t143;
        double t243 = t102+t142;

        double t247 = -t212;
        double t250 = t8*t218;
        double t256 = Itzz*Iwyy*t126*t142;

        double t258 = t246*4.0;
        double t259 = t8*t215;

        double t263 = t68*t177*dth*dps;
        double t267 = rw*t5*t210;
        double t274 = t106+t186;
        double t276 = t5*t246;

        double t279 = t3*t236*4.0;
        double t280 = g*t3*t10*t190*4.0;
        double t282 = rw*t8*t10*t190*4.0;
        double t283 = Iwyy*t87*t126*t145;
        double t286 = t3*t5*t70*t172;
        double t290 = t11*t124*t186;
        double t291 = t57*t61*t172;
        double t292 = t70*t124*t128;

        double t299 = t61*t70*t172*dth;
        double t301 = t288*2.0;
        double t305 = t4*t9*t125*t169;
        double t306 = t288*dph;
        double t307 = t275/2.0;
        double t311 = t16+t149+t153;
        double t312 = t10*t70*t96*t174*dth;
        double t315 = t10*t67*t70*t172;
        double t330 = t6*t11*t100*t124*t126;
        double t331 = t6*t11*t102*t124*t126;
        double t338 = t10*t54*t70*t174*dph;
        double t339 = t51*t287;
        double t346 = t61*t63*t75*t175;
        double t350 = rw*t30*t59*t62*t190;
        double t358 = t215*t236;
        double t368 = t104+t159+t223;
        double t370 = Iwyy+t215+t218;
        double t380 = Itxx+Ityy+t16+t81+t88+t130+t166;
        double t220 = Iwyy*t13*t185;
        double t221 = Iwyy*mt*rf*t185;
        double t222 = Iwyy*t14*t185;

        double t238 = rw*t228;
        double t248 = -t229;
        double t249 = -t231;
        double t251 = t234*2.0;
        double t252 = t235*2.0;
        double t255 = rw*t227*4.0;
        double t264 = -t234;
        double t266 = rw*t5*t209;
        double t268 = rw*t5*t211;
        double t269 = rw*t10*t209;
        double t270 = rw*t211*dth*dom;
        double t271 = Itzz+t247;
        double t272 = t100+t188;
        double t273 = t105+t185;

        double t281 = t235/4.0;
        double t284 = t28+t229;
        double t285 = t33+t230;
        double t289 = -t280;

        double t297 = t54*t239;
        double t298 = t56*t239;
        double t300 = t59*t237*dps;
        double t308 = t276/2.0;

        double t316 = rw*t10*t72*t211;
        double t317 = t3*t10*t241*4.0;
        double t318 = t126+t227;
        double t319 = -t307;
        double t321 = t126+t228;
        double t324 = t147+t203+1.0;
        double t327 = -t299;
        double t329 = rw*t8*t10*t51*t211;
        double t332 = t291/4.0;
        double t334 = -t305;
        double t337 = t53+t123+t202;
        double t343 = -t315;

        double t351 = -t338;
        double t353 = t74*t311;
        double t354 = t190*t242;
        double t356 = -t346;
        double t361 = Iwyy*t190*t243;
        double t366 = rw*t10*t190*t233;

        double t372 = rw*t5*t204*t243;
        double t374 = t228*t290;
        double t375 = rw*t2*t3*t190*t243;
        double t376 = t158+t187+t200;
        double t379 = rw*t7*t8*t190*t243;
        double t387 = t2*t126*t190*t243;
        double t398 = t56*t380;
        double t399 = Itxx+Ityy+t18+t130+t166+t247;
        double t401 = t279+t282;
        double t406 = Ifxx+t80+t87+t137+t138+t239;
        double t409 = t6*t11*t124*t228*t243;
        double t426 = t30*t59*t61*t370;
        double t428 = t51*t124*t244*t274;
        double t438 = t212+t380;
        double t458 = t17+t22+t25+t82+t89+t131+t167+t229;
        double t495 = t160+t161+t162+t195+t196+t197;
        double t663 = t118+t178+t179+t180+t181+t182+t183+t208+t217+t232+t240;
        double t260 = Iwyy+t238;
        double t265 = -t251;
        double t295 = t28+t248;
        double t296 = t34+t249;
        double t303 = t268/4.0;
        double t309 = t285*dps;
        double t314 = t5*t271*dth;

        double t325 = -t297;
        double t326 = -t298;
        double t328 = t67*t268*dph;
        double t335 = t3*t10*t271;
        double t340 = t5*t8*t284;
        double t345 = Iwyy*t8*t10*t87*t238;
        double t359 = g*t8*t318*4.0;
        double t360 = t70*t324;
        double t362 = Iwyy*t2*t3*t321;
        double t363 = t190*t273;
        double t364 = t5*t128*t284;
        double t365 = t3*t184*t268*dph;
        double t369 = t8*t361;
        double t371 = t214+t269;
        double t373 = -t366;
        double t378 = Iwyy+t213+t266;
        double t381 = rw*t59*t71*t321;
        double t382 = rw*t2*t3*t7*t8*t321;
        double t383 = Itzz*t59*t100*t321;
        double t384 = Itzz*t59*t102*t321;
        double t386 = rw*t5*t204*t272;
        double t388 = -t375;
        double t389 = rw*t2*t8*t190*t272;
        double t390 = rw*t3*t7*t190*t272;
        double t393 = t3*t58*t126*t321;
        double t394 = t3*t71*t126*t321;
        double t395 = t70*t124*t337;
        double t396 = t7*t126*t190*t272;
        double t404 = t199+t354;

        double t408 = t399*dom;
        double t410 = t19+t34+t150+t152+t255;
        double t411 = t398/4.0;
        double t412 = t3*t6*t11*t100*t124*t321;
        double t413 = rw*t3*t100*t190*t321;
        double t414 = t3*t6*t11*t102*t124*t321;
        double t415 = t8*t399;
        double t416 = rw*t3*t102*t190*t321;
        double t417 = t53*t401;
        double t419 = t19+t22+t25+t131+t167+t248;
        double t420 = t6*t11*t124*t228*t272;
        double t421 = t59*t101*t190*t321;

        double t423 = t59*t105*t190*t321;
        double t424 = (Iwyy*t399)/2.0;
        double t433 = t69*t406;

        double t453 = t3*t10*t438;
        double t454 = t8*t190*t243*t321;
        double t455 = t69*t438;
        double t461 = (t2*t126*t399)/2.0;
        double t462 = (t7*t126*t399)/2.0;
        double t468 = t3*t7*t190*t243*t321;
        double t470 = t8*t190*t272*t321;
        double t478 = t2*t3*t190*t272*t321;
        double t480 = t267+t399;
        double t482 = t69*t458*dth*dom;
        double t488 = t3*t56*t458*dph*dth;
        double t494 = t3*t56*t184*t438*dph;

        double t509 = (t2*t3*t321*t399)/2.0;
        double t511 = (t3*t7*t321*t399)/2.0;
        double t514 = t3*t56*t458*(t96-dom);
        double t516 = t163+t164+t165+t220+t221+t222;
        double t533 = t236*t495;
        double t563 = Itxx+Ityy+t19+t34+t77+t81+t88+t90+t94+t116+t132+t205+t235+t264;
        double t569 = t190*t272*t495;
        double t606 = (t399*t495)/2.0;
        double t680 = Itxx+Ityy+t16+t31+t35+t130+t148+t166+t235+t268+t291+t398;
        double t691 = t236*t663;
        double t717 = t190*t243*t663;
        double t720 = t2*t3*t321*t663;
        double t721 = t3*t7*t321*t663;
        double t723 = t190*t272*t663;
        double t302 = t8*t260;
        double t323 = t296*dps;
        double t341 = t335*dph;
        double t342 = -t314;
        double t344 = t5*t295*dom;
        double t347 = Iwyy*t72*t260;
        double t352 = t5*t8*t295*2.0;

        double t385 = t5*t11*t103*t124*t260;
        double t391 = t3*t371;

        double t397 = t3*t378*dph*dth;
        double t403 = t246+t335;
        double t405 = -t396;

        double t429 = -t414;
        double t430 = Itzz*t10*t404;
        double t431 = t419*dom;
        double t432 = -t424;

        double t435 = t288+t340;
        double t436 = t415/2.0;

        double t440 = t72*t410;
        double t445 = t3*t419*dph*dth;
        double t448 = t2*t126*t404;
        double t449 = -Itzz*t10*(t198-t363);
        double t450 = Iwyy*t10*t113*t415;
        double t456 = t236+t390;
        double t459 = t59*t433*2.0;

        double t465 = -t461;
        double t466 = -t462;
        double t471 = -t7*t126*(t198-t363);
        double t472 = -t124*dph*(t226-t360);
        double t473 = t236*t404;
        double t474 = t270+t395;
        double t475 = t238*t404;

        double t485 = t214+t433;
        double t489 = t236*(t198-t363);

        double t492 = (t8*t480)/2.0;
        double t493 = t309+t408;
        double t498 = t331+t413;
        double t499 = t3*t7*t321*t404;
        double t500 = t237+t455;

        double t506 = t190*t272*t404;

        double t513 = t455*(t52-t169);
        double t515 = -t190*t243*(t198-t363);

        double t518 = t393+t394;
        double t520 = -t5*(t246-t453);
        double t521 = t306+t327+t364;
        double t522 = Iwyy*t3*t5*t113*t321*t415;

        double t527 = (t7*t126*(t275-t415))/2.0;

        double t529 = rw*t2*t8*t516;
        double t540 = t151+t154+t225+t326+t353;
        double t543 = t236*(t275-t415)*(-1.0/2.0);
        double t546 = (t399*t404)/2.0;
        double t551 = t236*t516;

        double t573 = t3*t7*t321*(t275-t415)*(-1.0/2.0);

        double t577 = t421+t454;
        double t581 = t190*t243*t516;
        double t592 = t22+t25+t78+t79+t82+t89+t91+t92+t95+t134+t135+t219+t252+t265;
        double t600 = (t51*t67*t563)/2.0;

        double t619 = t468+t478;
        double t621 = t409+t509;
        double t622 = ((t198-t363)*(t275-t415))/2.0;
        double t634 = (t399*t516)/2.0;

        double t657 = t238*(t423-t470);

        double t671 = -rw*t3*t7*(t420-t511);
        double t676 = rw*t2*t8*(t420-t511);

        double t694 = ((t275-t415)*(t330-t416))/2.0;

        double t696 = (Itzz*t3*t5*t680)/4.0;
        double t703 = (t7*t126*t680)/4.0;
        double t712 = (Iwyy*t199*t680)/4.0;
        double t713 = (t215*t680)/4.0;

        double t725 = (Iwyy*t8*t238*t680)/4.0;
        double t741 = (t190*t243*t680)/4.0;
        double t744 = (t190*t272*t680)/4.0;
        double t759 = (t399*t680)/8.0;
        double t763 = (t404*t680)/4.0;

        double t792 = t663*(t275-t415)*(-1.0/2.0);

        double t810 = t680*(t275-t415)*(-1.0/8.0);
        double t825 = t109+t110+t111+t112+t113+t114+t115+t171+t201+t281+t303+t332+t381+t411;

        double t400 = t259+t302;
        double t402 = t391/2.0;
        double t418 = t403*dph*dom;
        double t441 = t190*t243*t302;
        double t442 = t435*dph;
        double t443 = t2*t3*t302*t321;
        double t444 = t3*t7*t302*t321;
        double t446 = t190*t272*t302;

        double t464 = -t459;
        double t469 = Iwyy*t456;

        double t486 = t10*t474*2.0;
        double t487 = t302*t404;
        double t491 = t3*t485;

        double t501 = t321*t436;
        double t502 = t3*t493*dph*2.0;

        double t510 = t292+t341+t342;
        double t512 = t8*t500;
        double t526 = t387+t405;
        double t530 = t521*dne*2.0;
        double t531 = Iwyy*t8*t518;
        double t532 = rw*t2*t518;
        double t534 = t383+t430;

        double t539 = rw*t3*t7*t518;
        double t548 = t384+t449;
        double t549 = t290+t492;
        double t552 = t388+t456;
        double t553 = t372+t466;
        double t555 = t236*t518;
        double t557 = t386+t465;
        double t558 = t59*t540*2.0;
        double t559 = t3*t10*t190*t518;
        double t575 = t190*t272*t498;

        double t598 = t448+t471;
        double t599 = t319+t389+t436;
        double t604 = t96*t592;
        double t611 = t391+t520;
        double t620 = t238*t577;
        double t626 = t254+t334+t373+t397+Tomega;
        double t631 = Iwyy*t621;
        double t642 = (t399*t518)/2.0;

        double t652 = rw*t2*t3*t619;
        double t654 = rw*t2*t8*t619;
        double t655 = rw*t3*t7*t619;
        double t656 = rw*t7*t8*t619;

        double t700 = t506+t515;
        double t714 = -Iwyy*(t499+t2*t3*t321*(t198-t363));
        double t749 = t302*(t499+t2*t3*t321*(t198-t363));
        double t802 = t569+t581;
        double t822 = t356+t759;
        double t824 = t412+t741;
        double t827 = Iwyy*t825;
        double t828 = t429+t744;
        double t832 = -t2*t126*(t346-t759);

        double t836 = ((t499+t2*t3*t321*(t198-t363))*(t275-t415))/2.0;
        double t840 = t7*t126*(t346-t759);
        double t842 = -rw*t2*t8*(t414-t744);
        double t867 = t302*(t346-t759);
        double t898 = t551+t721;
        double t930 = t634+t723;

        double t467 = t236*t400;
        double t476 = t250+t400;
        double t497 = t491/2.0;
        double t519 = t510*dne*4.0;

        double t535 = -t530;
        double t536 = rw*t2*t526;
        double t538 = t8*t532;
        double t541 = Iwyy*t534;
        double t542 = dth*(t344-t442)*(-1.0/2.0);
        double t545 = rw*t3*t7*t526;

        double t550 = Iwyy*t548;
        double t556 = Itzz*t10*t549;
        double t560 = rw*t7*t553;
        double t561 = t286+t512;

        double t565 = rw*t2*(t276-t491)*(-1.0/2.0);
        double t566 = t2*t126*t549;
        double t567 = t7*t126*t549;
        double t568 = rw*t2*t557;

        double t571 = t238*t534;
        double t574 = t287+t316+t464;
        double t578 = t3*t5*t87*t552;

        double t583 = t238*t548;
        double t584 = t236*t549;
        double t589 = t8*t321*t526;
        double t590 = t374+t501;
        double t593 = t236*(t276-t491)*(-1.0/2.0);
        double t601 = t302*t549;
        double t603 = Iwyy*t599;
        double t607 = t362+t532;
        double t608 = rw*t2*t598;
        double t616 = rw*t3*t7*t598;
        double t618 = t190*t243*t549;
        double t623 = t2*t3*t321*t549;
        double t624 = t3*t7*t321*t549;
        double t627 = t190*t272*t549;
        double t630 = -t620;
        double t633 = t190*t243*(t276-t491)*(-1.0/2.0);

        double t636 = t3*t7*t321*(t276-t491)*(-1.0/2.0);
        double t637 = t361*(t276-t491)*(-1.0/2.0);
        double t643 = (t2*t3*t321*(t276-t491))/2.0;

        double t649 = t190*t272*(t276-t491)*(-1.0/2.0);
        double t650 = t8*t631;

        double t658 = t236*t598;
        double t662 = t3*t5*t113*t611;
        double t666 = (t2*t126*t611)/2.0;
        double t667 = (t7*t126*t611)/2.0;
        double t668 = (t399*t534)/2.0;
        double t672 = Itzz*t3*t74*t190*t598;
        double t678 = t399*(t276-t491)*(-1.0/4.0);
        double t679 = (t399*t548)/2.0;
        double t681 = (t236*t611)/2.0;
        double t682 = (t238*t611)/2.0;
        double t685 = (t3*t10*t190*t611)/2.0;
        double t689 = -rw*t2*(t446+(t2*t126*(t275-t415))/2.0);

        double t693 = t379+t599;
        double t698 = (t190*t243*t611)/2.0;
        double t708 = (t190*t272*t611)/2.0;

        double t716 = rw*t2*t700;
        double t719 = (t399*t598)/2.0;
        double t728 = ((t275-t415)*(t276-t491))/4.0;

        double t731 = (t399*t611)/4.0;
        double t733 = (t404*t611)/2.0;
        double t735 = -Iwyy*(t382-t402+(t5*(t246-t453))/2.0);

        double t747 = t526*(t276-t491)*(-1.0/2.0);

        double t756 = (t400*t680)/4.0;
        double t760 = (t598*(t275-t415))/2.0;
        double t767 = t323+t431+t604;
        double t813 = (t549*t611)/2.0;

        double t815 = rw*t2*t8*t802;
        double t816 = rw*t3*t7*t802;
        double t817 = rw*t7*t8*t802;
        double t820 = (t611*(t276-t491))/4.0;
        double t830 = (t526*t680)/4.0;
        double t834 = t117+t133+t206+t325+t343+t440+t558;

        double t843 = ((t441+t527)*(t276-t491))/2.0;
        double t849 = t238*t824;
        double t850 = (t549*t680)/4.0;
        double t852 = ((t276-t491)*(t446+(t2*t126*(t275-t415))/2.0))/2.0;
        double t864 = t539+t703;
        double t870 = t531+t714;
        double t896 = (t611*t663)/2.0;
        double t905 = Itzz*t3*t5*(t696-t10*t113*(t276-t491));
        double t916 = rw*t2*t8*t898;
        double t917 = rw*t3*t7*t898;

        double t929 = t680*(t446+(t2*t126*(t275-t415))/2.0)*(-1.0/4.0);
        double t935 = rw*t2*t8*t930;
        double t936 = rw*t3*t7*t930;
        double t971 = -Iwyy*t8*(-t575+t642+t190*t243*(t330-t416));

        double t976 = -Itzz*t3*t5*(-t575+t642+t190*t243*(t330-t416));

        double t979 = -rw*t3*t7*(-t575+t642+t190*t243*(t330-t416));
        double t980 = rw*t2*(-t575+t642+t190*t243*(t330-t416));
        double t982 = t655+t824;

        double t989 = t652+t828;
        double t990 = t671+t822;

        double t1086 = -rw*t2*(t832+t238*(t414-t744));
        double t1090 = Iwyy*t8*(t832+t238*(t414-t744));
        double t1092 = rw*t2*(t832+t238*(t414-t744));

        double t544 = t8*t536;
        double t579 = t52*t561;
        double t582 = t574*dom;

        double t588 = t283+t541;

        double t594 = Itzz*t10*t399*t476*(-1.0/2.0);
        double t595 = t191+t565;
        double t596 = t256+t550;
        double t597 = Iwyy*t590;
        double t609 = t8*t603;
        double t610 = t361+t536;

        double t613 = t2*t126*t590;
        double t614 = t7*t126*t590;
        double t615 = t8*t608;
        double t628 = -t618;
        double t629 = t336+t556;

        double t639 = rw*t2*t3*t607;
        double t640 = rw*t7*t8*t607;
        double t648 = rw*t7*(t331-t545);

        double t665 = t302*t590;
        double t673 = -t667;

        double t701 = Itzz*t10*t693;

        double t710 = -Iwyy*(t308-t497+rw*t3*t7*(t198-t363));
        double t722 = t8*t716;

        double t736 = -rw*t7*(t475-t567);
        double t742 = rw*t2*(t566-t238*(t198-t363));
        double t743 = rw*t7*t8*t87*(t475-t567);
        double t748 = t385+t682;
        double t751 = t473+t623;
        double t754 = t473+t633;
        double t757 = t444+t666;
        double t758 = t3*t74*t87*t190*(t475-t567);
        double t761 = t495+t608;

        double t765 = Itzz*t3*t74*t190*(t566-t238*(t198-t363));

        double t774 = t3*t767;
        double t775 = t489+t649;
        double t778 = t190*t243*(t566-t238*(t198-t363));

        double t797 = t236*(t489-t624);
        double t799 = t3*t10*t190*(t489-t624);
        double t800 = t559+t589;

        double t835 = t253+t333+t418+t428+t542+Tneta;
        double t841 = rw*t7*(t616-(t7*t126*(t276-t491))/2.0);
        double t846 = -t843;
        double t847 = (Iwyy*t834)/4.0;
        double t853 = (Itzz*t10*t834)/4.0;
        double t854 = ((t275-t415)*(t475-t567))/2.0;

        double t859 = (t2*t126*t834)/4.0;
        double t860 = (t7*t126*t834)/4.0;

        double t865 = ((t566-t238*(t198-t363))*(t275-t415))/2.0;
        double t873 = t8*t321*(t627-(t399*(t198-t363))/2.0);
        double t875 = rw*t7*t864;
        double t876 = t584+t678;
        double t879 = (t236*t834)/4.0;
        double t880 = (t238*t834)/4.0;
        double t887 = t432+t560+t568;
        double t888 = rw*t2*t8*t870;
        double t889 = rw*t3*t7*t870;

        double t903 = t238*(t685-(t8*t321*(t275-t415))/2.0);
        double t904 = t573+t708;
        double t906 = ((t275-t415)*(t489-t624))/2.0;
        double t907 = (t190*t243*t834)/4.0;
        double t910 = (t2*t3*t321*t834)/4.0;
        double t911 = (t3*t7*t321*t834)/4.0;
        double t913 = (t190*t272*t834)/4.0;

        double t942 = t658+t747;
        double t947 = t643+t763;
        double t952 = t555+t830;
        double t959 = -t7*t126*(t636+(t680*(t198-t363))/4.0);
        double t985 = t593+t850;
        double t993 = t681+t810;

        double t995 = rw*t7*t8*t982;
        double t996 = Iwyy*t990;
        double t997 = t97+t289+t312+t339+t445+t488+t513+t519;
        double t1000 = t3*t5*t87*t982;
        double t1004 = Itzz*t3*t5*t989;
        double t1011 = t399*(t636+(t680*(t198-t363))/4.0)*(-1.0/2.0);
        double t1015 = (t680*t834)/1.6e+1;
        double t1027 = -Itzz*t10*(-t654+t698+(t2*t3*t321*(t275-t415))/2.0);
        double t1030 = -rw*t2*t3*(-t654+t698+(t2*t3*t321*(t275-t415))/2.0);
        double t1031 = t543+t676+t731;
        double t1033 = ((t636+(t680*(t198-t363))/4.0)*(t275-t415))/2.0;
        double t1073 = t840+t849;
        double t1074 = t631+t980;
        double t1091 = t8*t1086;
        double t1104 = t263+t359+t365+t482+t486+t494+t502+t535+t600;
        double t587 = -t582;
        double t602 = t8*t597;

        double t625 = t2*t126*t595;
        double t638 = Iwyy*t629;
        double t644 = rw*t2*t3*t610;
        double t646 = rw*t7*t8*t610;

        double t661 = t10*t87*t629;
        double t664 = t7*t126*t629;
        double t669 = t2*t126*t629;
        double t697 = t190*t243*t629;
        double t706 = t190*t272*t629;
        double t711 = (t399*t588)/2.0;
        double t718 = (t399*t596)/2.0;
        double t745 = Itzz*t8*t742;
        double t766 = rw*t2*t757;
        double t768 = t2*t126*t751;

        double t771 = rw*t3*t7*t757;
        double t772 = rw*t7*t8*t761;
        double t776 = rw*t2*t3*t761;
        double t777 = Iwyy*t8*(t443+t673);
        double t781 = t236*t751;

        double t783 = t3*t10*t190*t751;
        double t784 = t238*t754;
        double t785 = rw*t2*t8*t775;
        double t788 = t236*t757;
        double t793 = t302*t751;
        double t798 = t190*t243*t748;

        double t804 = t190*t272*t748;
        double t806 = t238*t775;
        double t807 = t546+t628;
        double t809 = t190*t272*t751;
        double t812 = t190*t243*t757;
        double t818 = Itzz*t3*t5*t800;

        double t829 = t441+t527+t544;
        double t845 = (t399*t757)/2.0;

        double t866 = -t860;
        double t871 = t614+t630;
        double t881 = -t880;
        double t884 = t613+t657;

        double t899 = t2*t126*t876;
        double t900 = t7*t126*t876;
        double t902 = t578+t701;
        double t908 = -t903;

        double t912 = t3*t5*t10*t30*t887;
        double t914 = -t910;
        double t920 = rw*t3*t7*t904;
        double t922 = t302*t876;
        double t923 = t443+t538+t673;
        double t938 = t190*t243*(t347-t847);
        double t945 = rw*t2*t942;

        double t949 = Iwyy*t947;
        double t950 = rw*t2*(t859-t302*(t198-t363));
        double t953 = t2*t126*t947;
        double t955 = t328+t351+t514+t774;
        double t958 = t236*(t859-t302*(t198-t363));

        double t961 = rw*t2*t8*t952;
        double t962 = rw*t7*t8*t952;
        double t963 = t238*t947;
        double t965 = t3*t5*t87*t952;
        double t970 = t302*t947;
        double t975 = t190*t243*(t859-t302*(t198-t363));
        double t984 = t190*t272*t947;
        double t987 = Itzz*t985;
        double t988 = Iwyy*t985;
        double t1001 = t8*t996;
        double t1002 = t2*t126*t985;
        double t1003 = t7*t126*t985;
        double t1007 = (t399*(t859-t302*(t198-t363)))/2.0;
        double t1008 = (t399*t947)/2.0;
        double t1013 = t662+t853;
        double t1014 = t238*t993;

        double t1018 = t190*t243*t985;
        double t1019 = t190*t272*t985;
        double t1023 = t656+t904;
        double t1024 = t622+t913;
        double t1036 = Iwyy*t1031;
        double t1042 = t728+t879;
        double t1043 = t985*(t275-t415)*(-1.0/2.0);
        double t1058 = -t7*t126*(t911+(t611*(t198-t363))/2.0);
        double t1064 = -t236*(t911+(t611*(t198-t363))/2.0);
        double t1066 = t799+t873;
        double t1071 = t813+t879;
        double t1076 = t663+t736+t742;
        double t1078 = Iwyy*t8*t1073;
        double t1083 = rw*t2*t3*t1074;
        double t1084 = rw*t7*t8*t1074;
        double t1095 = t399*(t911+(t611*(t198-t363))/2.0)*(-1.0/2.0);
        double t1101 = -rw*t7*(-t722+t907+(t404*(t275-t415))/2.0);
        double t1106 = -rw*t2*(-t719+t778+t190*t272*(t475-t567));
        double t1120 = t820+t1015;
        double t1170 = t1000+t1027;
        double t1173 = t639+t713+t827+t875;
        double t1182 = t979+t1073;

        double t1204 = t10*t87*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567)));
        double t1206 = -rw*t7*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567)));
        double t659 = -t644;
        double t674 = -t664;
        double t675 = t37+t300+t587;
        double t702 = t345+t638;
        double t707 = -t697;
        double t838 = rw*t7*t829;
        double t848 = -t845;
        double t862 = t8*t321*t807;
        double t863 = t583+t669;
        double t882 = Iwyy*t8*t871;
        double t883 = rw*t7*t871;
        double t892 = Itzz*t3*t5*t871;
        double t894 = Iwyy*t8*t884;
        double t895 = rw*t2*t884;
        double t901 = Itzz*t3*t5*t884;

        double t924 = t238*t902;
        double t925 = rw*t7*t923;
        double t933 = t679+t706;
        double t943 = t487+t866;
        double t948 = t8*t945;
        double t957 = t955*dth;
        double t991 = t672+t818;

        double t1006 = -t1003;
        double t1009 = -t1008;
        double t1010 = t601+t881;
        double t1016 = Itzz*t10*t1013;
        double t1020 = t756+t771;
        double t1028 = rw*t2*t1024;
        double t1029 = Itzz*t10*t1023;
        double t1034 = rw*t3*t7*t1024;
        double t1035 = t665+t908;
        double t1040 = t8*t1036;
        double t1044 = t238*t1042;
        double t1046 = t733+t914;
        double t1067 = Iwyy*t1066;
        double t1068 = rw*t2*t1066;

        double t1080 = t2*t126*t1071;
        double t1081 = t7*t126*t1071;
        double t1082 = -rw*t2*(t806-t899);
        double t1093 = t236*t1071;
        double t1094 = t661+t987;
        double t1096 = t190*t243*t1071;
        double t1097 = t190*t272*t1071;
        double t1107 = t8*t1106;

        double t1110 = t852+t958;
        double t1111 = -rw*t2*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624));
        double t1113 = -rw*t3*t7*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624));
        double t1117 = -t3*t10*t190*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624));
        double t1125 = t87*t1120;
        double t1126 = -Itzz*t10*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673));
        double t1127 = -Iwyy*t8*(t798+(t399*(t443+t673))/2.0+(t498*(t275-t415))/2.0);
        double t1130 = -rw*t2*t3*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673));
        double t1134 = t2*t126*t1120;
        double t1135 = t7*t126*t1120;
        double t1141 = t190*t243*t1120;
        double t1142 = t190*t272*t1120;
        double t1146 = ((t275-t415)*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624)))/2.0;
        double t1148 = (t399*t1120)/2.0;
        double t1152 = -rw*t3*t7*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624));
        double t1156 = t953+t959;
        double t1159 = ((t276-t491)*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673)))/2.0;

        double t1163 = Itzz*(t1002-t238*(t636+(t680*(t198-t363))/4.0));
        double t1165 = rw*t2*(t1002-t238*(t636+(t680*(t198-t363))/4.0));
        double t1166 = t625+t710+t776+t841;
        double t1168 = t533+t637+t816+t945;
        double t1180 = rw*t7*t8*(t650+rw*t2*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)));
        double t1183 = rw*t7*t1182;
        double t1196 = -rw*t7*t8*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0));
        double t1211 = -rw*t7*(t798+(t399*(t443+t673))/2.0+t8*t980+(t498*(t275-t415))/2.0);
        double t1213 = -rw*t2*t3*(t798+(t399*(t443+t673))/2.0+t8*t980+(t498*(t275-t415))/2.0);
        double t1219 = t797+t1011+t1019;
        double t1233 = -rw*t2*t3*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)));
        double t1244 = -rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)));
        double t1291 = -Iwyy*(t1120+rw*t2*t8*(t636+(t680*(t198-t363))/4.0)+rw*t3*t7*(t911+(t611*(t198-t363))/2.0));
        double t1300 = -rw*t7*(t961+t236*(t443+t673)+(t680*(t441+t527))/4.0+rw*t3*t7*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673)));
        double t1323 = t842+t920+t993+t995+t1030;

        double t686 = t675*dph*2.0;
        double t737 = t3*t10*t190*t702;
        double t750 = t190*t243*t702;
        double t752 = t190*t272*t702;
        double t844 = -t838;
        double t855 = t571+t674;
        double t872 = rw*t2*t863;
        double t886 = t10*t87*t863;
        double t928 = -t925;
        double t931 = t668+t707;
        double t940 = rw*t2*t933;
        double t954 = t236*t943;
        double t973 = t190*t272*t943;
        double t981 = t358+t469+t648+t659;
        double t999 = (t399*t943)/2.0;
        double t1012 = Itzz*t1010;
        double t1022 = Iwyy*t8*t1020;
        double t1025 = t594+t924;
        double t1032 = t190*t243*t1010;
        double t1038 = t190*t272*t1010;
        double t1045 = t615+t943;
        double t1047 = -t1044;
        double t1048 = Iwyy*t1046;
        double t1052 = t758+t892;
        double t1053 = t2*t126*t1046;

        double t1057 = t765+t901;
        double t1060 = t236*t1046;
        double t1061 = t238*t1046;
        double t1069 = t190*t272*t1046;
        double t1070 = -t1068;
        double t1087 = -t1081;
        double t1088 = (t399*t1046)/2.0;
        double t1089 = t8*t1082;
        double t1099 = t302*t1094;
        double t1112 = t8*t1111;
        double t1114 = rw*t2*t1110;
        double t1116 = Itzz*t1113;
        double t1123 = t694+t804+t848;
        double t1138 = -t1135;
        double t1139 = t597+t883+t895;
        double t1149 = -t1148;

        double t1157 = rw*t2*t1156;
        double t1160 = t963+t1006;
        double t1167 = Itzz*t3*t5*t1166;
        double t1171 = rw*t7*t8*t1168;
        double t1172 = t1004+t1029;
        double t1184 = t533+t720+t1111;
        double t1185 = t8*t1183;
        double t1186 = -t1183;
        double t1208 = t965+t1126;
        double t1210 = t781+t1009+t1018;
        double t1216 = rw*t2*(t1134-t302*(t636+(t680*(t198-t363))/4.0));
        double t1221 = Iwyy*t1219;
        double t1222 = rw*t2*t1219;
        double t1225 = t7*t126*t1219;
        double t1227 = t302*t1219;
        double t1240 = t785+t1034+t1042;
        double t1246 = t8*t1244;
        double t1269 = t190*t243*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624));

        double t1271 = -rw*t7*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567)));
        double t1272 = t905+t1016+t1125;
        double t1275 = t906+t1095+t1097;
        double t1294 = -rw*t3*t7*(Itzz*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624))+t3*t5*t87*t863);
        double t1297 = t788+t929+t962+t1130;
        double t1303 = t1033+t1064+t1142;
        double t1315 = -t7*t126*(-t1093+t1148+(t985*(t275-t415))/2.0);
        double t1316 = t7*t126*(-t1093+t1148+(t985*(t275-t415))/2.0);
        double t1324 = Iwyy*t1323;

        double t868 = Itzz*t10*t855;
        double t869 = rw*t7*t855;
        double t877 = Itzz*t3*t5*t855;
        double t932 = rw*t7*t931;
        double t941 = -t940;
        double t998 = t3*t5*t87*t981;
        double t1039 = -t1032;
        double t1041 = t3*t5*t87*t1025;
        double t1049 = rw*t7*t1045;
        double t1055 = rw*t7*t1052;
        double t1063 = rw*t2*t1057;
        double t1115 = Itzz*t1112;
        double t1128 = rw*t2*t1123;
        double t1129 = Iwyy*t8*t1123;
        double t1137 = rw*t3*t7*t1123;
        double t1140 = t603+t646+t689+t844;
        double t1158 = t8*t1157;

        double t1164 = t87*t1160;
        double t1169 = t350+t743+t745+t1012;
        double t1175 = t640+t735+t766+t928;
        double t1188 = rw*t2*t3*t1184;
        double t1189 = rw*t7*t8*t1184;
        double t1194 = t760+t973+t975;
        double t1212 = Iwyy*t1210;
        double t1215 = t886+t1163;
        double t1217 = t2*t126*t1210;
        double t1218 = -t1216;
        double t1223 = t8*t1222;
        double t1224 = t302*t1210;
        double t1237 = t865+t1007+t1038;
        double t1241 = Iwyy*t1240;
        double t1247 = t749+t1053+t1058;
        double t1253 = t793+t1061+t1087;
        double t1268 = rw*t2*t3*(t836-t1069+t190*t243*(t911+(t611*(t198-t363))/2.0));

        double t1276 = Iwyy*t1275;
        double t1277 = t238*t1272;
        double t1278 = rw*t2*t1275;
        double t1279 = rw*t3*t7*t1275;
        double t1282 = t718+t752+t1206;
        double t1285 = -rw*t7*(t1160+rw*t3*t7*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624)));
        double t1296 = rw*t7*t8*(t712+t889-t949-t1157);
        double t1299 = rw*t2*t1297;
        double t1304 = -t238*(-t1060+t1141+(t947*(t275-t415))/2.0);
        double t1305 = t238*t1303;
        double t1311 = t1043+t1093+t1149;
        double t1321 = t1152+t1210;
        double t1336 = -rw*t7*(-t1088+t1096+(t751*(t275-t415))/2.0+rw*t2*t8*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)));
        double t1337 = t894+t1067+t1271;
        double t1344 = t996+t1083+t1092+t1186;
        double t874 = -t869;
        double t1133 = -t1128;
        double t1144 = Itzz*t10*t1140;
        double t1154 = t8*t321*t1140;
        double t1155 = t450+t932+t941;
        double t1181 = t3*t10*t190*t1175;

        double t1197 = rw*t2*t1194;
        double t1198 = rw*t7*t1194;
        double t1220 = rw*t2*t8*t1215;
        double t1236 = t854+t999+t1039;
        double t1238 = rw*t2*t1237;
        double t1239 = rw*t3*t7*t1237;
        double t1242 = (t680*t1194)/4.0;
        double t1248 = rw*t2*t1247;
        double t1249 = rw*t3*t7*t1247;
        double t1250 = t236*t1247;
        double t1252 = (t399*t1247)/2.0;
        double t1254 = Itzz*t1253;
        double t1258 = t522+t737+t1055+t1063;
        double t1262 = t190*t272*t1253;
        double t1280 = -t1278;
        double t1287 = -Itzz*(t347+t529-t772-t847-t950+t1049);
        double t1306 = t868+t1116+t1164;
        double t1322 = rw*t7*t1321;
        double t1325 = t867+t1014+t1091+t1137;
        double t1327 = t1112+t1253;
        double t1345 = Itzz*t3*t5*t1344;
        double t1352 = -Itzz*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)));
        double t1356 = rw*t7*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)));
        double t1362 = t1196+t1268+t1303;
        double t1364 = t3*t10*t190*(-t725+t917+t988+t1165+t1188+t1285);
        double t1368 = t1223+t1279+t1311;
        double t1382 = t1299+t1300+t1324;
        double t1153 = t702+t872+t874;
        double t1199 = t3*t1198;
        double t1200 = -t1197;

        double t1230 = t998+t1144;
        double t1251 = -t1250;
        double t1260 = t3*t5*t87*t1258;
        double t1263 = -t1262;
        double t1288 = t10*t1287;
        double t1301 = t1154+t1181;
        double t1307 = rw*t7*t8*t1306;
        double t1317 = t1107+t1236;
        double t1326 = Iwyy*t8*t1325;
        double t1329 = rw*t7*t1327;
        double t1338 = t877+t1115+t1254;
        double t1339 = t777+t888+t1048+t1248;
        double t1346 = t1036+t1084+t1133+t1211;
        double t1348 = t970+t1138+t1158+t1249;
        double t1360 = t1185+t1213+t1325;
        double t1367 = t1204+t1352;
        double t1369 = Iwyy*t1368;
        double t1381 = t1040+t1180+t1280+t1336;
        double t1394 = t1090+t1221+t1233+t1356;
        double t1232 = Itzz*t3*t5*t1230;
        double t1313 = rw*t2*t3*(t815+t938+t1200-(t495*(t275-t415))/2.0);
        double t1318 = rw*t7*t1317;
        double t1319 = rw*t2*t3*t1317;

        double t1331 = t846+t948+t954+t1199;
        double t1340 = rw*t2*t3*t1338;
        double t1341 = rw*t2*t3*t1339;
        double t1342 = t1167+t1288;
        double t1347 = -t10*t87*t1346;
        double t1349 = rw*t7*t1348;
        double t1350 = t1159+t1242+t1251;
        double t1372 = t1146+t1252+t1263+t1269;
        double t1384 = -t3*t10*t190*(-t896+t916-t1189+t1329+t236*(t347-t847)+rw*t2*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)));

        double t1332 = rw*t7*t1331;
        double t1374 = rw*t2*t1372;
        double t1375 = rw*t7*t1372;
        double t1386 = t1345+t1347;
        double t1389 = t922+t1047+t1089+t1239+t1246+t1319;
        double t1405 = -Itzz*(t1384+t8*t321*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))));
        double t1406 = t1022+t1218+t1291+t1296+t1341+t1349;
        double t1333 = -t1332;
        double t1376 = t3*t1375;
        double t1390 = Itzz*t1389;
        double t1408 = t1260+t1405;
        double t1396 = Itzz*(t1114+t1171+t1241+t1313+t1333-Iwyy*t8*(t467+rw*t3*t7*(t446+(t2*t126*(t275-t415))/2.0)));
        double t1409 = rw*t8*t1408;
        double t1414 = -rw*t7*(-t1224+t1315+t1376+t238*(-t1060+t1141+(t947*(t275-t415))/2.0)+rw*t2*t8*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0))));
        double t1410 = t1409*4.0;

        double et1 = t1409;
        double et2 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et3 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double t1421 = -1.0/(et1+et2+et3);
        double et4 = t1410;
        double et5 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et6 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double t1422 = -1.0/(et4+et5+et6);

        double et7 = t1410;
        double et8 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et9 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et10 = t1409;
        double et11 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et12 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et13 = t1409;
        double et14 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et15 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et16 = t1409;
        double et17 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et18 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et19 = t1422*(Itzz*t1344+t30*t74*t887+Itzz*rw*t8*t1139)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352)))+t835*t1421*(t1345-Itzz*t10*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+rw*t8*t1258);
        double et20 = t1104*t1422*(t912+Itzz*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))-Itzz*rw*t3*t8*t10*t190*t1076)+(t997*(Itzz*(-t725+t917+t988+t1165+t1188+t1285)+t10*t87*t1153+Itzz*rw*t72*t321*t1076))/(et7+et8+et9);
        double et21 = (t626*(Itzz*(t1001-t1222-t1322+rw*t2*t3*(t650+rw*t2*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))-Itzz*t10*t1155+rw*t8*t87*(-t602+t1068+rw*t7*(t783-t862))))/(et10+et11+et12);
        double et22 = (rw*t368*(-Itzz*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+t10*t87*(-t711+t750+rw*t2*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567))))+Itzz*rw*t8*(t882+Iwyy*(t783-t862)+rw*t2*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))))))/(et13+et14+et15)+rw*t376*t1421*(Itzz*t1394+t10*t87*t1282+Itzz*rw*t8*t1337);
        double et23 = (rw*t3*t52*(Itzz*(t1364+t8*t321*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t10*t1258))/(et16+et17+et18);
        double et24 = t1410;
        double et25 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et26 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et27 = t1409;
        double et28 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et29 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et30 = t1409;
        double et31 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et32 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et33 = t1409;
        double et34 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et35 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et36 = t1409;
        double et37 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et38 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et39 = ((t912-Itzz*t1346)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et24+et25+et26)+(t835*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+t3*t5*t87*t1346))/(et27+et28+et29)+t1104*t1422*(Itzz*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))-t30*t59*t61*t887)+(t626*(Itzz*t1381+Itzz*t3*t5*t1155))/(et30+et31+et32);
        double et40 = t997*t1422*(Itzz*(-t896+t916-t1189+t1329+t236*(t347-t847)+rw*t2*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)))+t3*t5*t87*t1153)+(rw*t368*(Itzz*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))-t3*t5*t87*(-t711+t750+rw*t2*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567))))))/(et33+et34+et35);
        double et41 = (rw*t376*(Itzz*(t1129-t1276-t1375+rw*t7*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+t3*t5*t87*t1282))/(et36+et37+et38)+rw*t3*t52*t1408*t1421;
        double et42 = t1410;
        double et43 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et44 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et45 = t1409;
        double et46 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et47 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et48 = t1409;
        double et49 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et50 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et51 = t87*(-t1093+t1148+(t985*(t275-t415))/2.0)+Itzz*t10*((t399*t1013)/2.0+(t629*(t275-t415))/2.0)-rw*t8*(Itzz*((t590*t834)/4.0-t549*(t685-(t8*t321*(t275-t415))/2.0))+Itzz*t8*t1070-Itzz*t3*t5*(t3*t10*t190*t629+t3*t5*t113*t321*t415)+rw*t7*t8*t87*(t783-t862))-rw*t2*t3*(Itzz*(-t1088+t1096+(t751*(t275-t415))/2.0)-Itzz*t3*t5*t931+Itzz*rw*t2*t8*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))+rw*t3*t7*(Itzz*t1275-t3*t5*t87*t933);
        double et52 = Itzz*t3*t5*(t236*t629+(t399*(t696-t10*t113*(t276-t491)))/2.0)+rw*t7*t8*(Itzz*t1152+Itzz*t1210+Itzz*t10*t931)+rw*t2*t8*(Itzz*t1219+t10*t87*t933);
        double et53 = t997*t1422*(t1099+t1220+t1277+t1294+t1307+t1340+rw*t72*t321*t1169)+((t87*t1360+Itzz*t10*t1025+rw*t8*(Itzz*t1035+Itzz*t8*t883+Itzz*t8*t895-rw*t30*t59*t61*t74*t204))*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et42+et43+et44)+t1104*t1422*(t1041+t1390+rw*t3*t8*t10*t190*t1169)+(t835*(t10*t1390+rw*t8*(t8*t1055+t8*t1063+Itzz*t3*t5*t1035+t3*t74*t190*t1012)+t3*t5*t87*t1360))/(et45+et46+et47)+(t626*(et51+et52))/(et48+et49+et50);
        double et54 = rw*t376*t1421*(t87*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+Itzz*t10*(Itzz*t10*t1237+t3*t5*t87*t1123)+rw*t8*(t87*(t8*t321*t1237+t3*t10*t190*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)))+Itzz*t3*t5*t1057+Itzz*rw*t7*t8*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))))+rw*t7*t8*t1367+rw*t2*t3*(Itzz*t1372-Itzz*t3*t5*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567))))-Itzz*t3*t5*(Itzz*t10*(t806-t899)+t3*t5*t87*(t832+t238*(t414-t744))));
        double et55 = rw*t368*t1421*(Itzz*(t1224+t1304+t1316)-t10*t87*(Itzz*t10*t1236+Itzz*t3*t5*(t798+(t399*(t443+t673))/2.0+(t498*(t275-t415))/2.0))+rw*t8*(-Itzz*(t8*t321*t1236-t3*t10*t190*t1253)+t3*t5*t87*t1052+Itzz*rw*t2*t8*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))))-t3*t5*t87*(Itzz*t10*(t784-t900)-t3*t5*t87*t1073)+rw*t2*t8*t1367-rw*t3*t7*(Itzz*t1372-Itzz*t3*t5*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567)))))+rw*t3*t122*t1421*(t8*t321*(t1041+t1390)-t3*t10*t190*(t1099+t1220+t1277+t1294+t1307+t1340));
        double et56 = t1410;
        double et57 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et58 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et59 = -Itzz*(-t8*t1324+rw*t2*t1362+rw*t7*(-t1060+t1141+(t947*(t275-t415))/2.0+rw*t2*t8*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0))+rw*t3*t7*(t836-t1069+t190*t243*(t911+(t611*(t198-t363))/2.0))));
        double et60 = rw*t8*(Itzz*(t8*t321*(t609-t1028+t1101+rw*t7*t8*(t369-t716))-t3*t10*t190*(Iwyy*t8*(t382-t402+(t5*(t246-t453))/2.0)-rw*t2*(t911+(t611*(t198-t363))/2.0)+rw*t7*(-t733+t910+rw*t2*t8*(t499+t2*t3*t321*(t198-t363)))-rw*t7*t8*(t8*t362+rw*t2*(t499+t2*t3*t321*(t198-t363)))))-t3*t5*t87*(rw*t7*(t3*t10*t190*t534+Itzz*t3*t8*t228*t243*t321)-rw*t2*(t3*t10*t190*t548+t3*t8*t87*t228*t272*t321)+Itzz*Iwyy*t3*t8*t74*t190));
        double et61 = t10*t87*(Itzz*t10*(t609-t1028+t1101+rw*t7*t8*(t369-t716))-Itzz*t3*t5*(t8*t469+rw*t7*(t754+rw*t3*t7*t700)+rw*t2*t775-rw*t2*t3*(t369-t716)))+Itzz*t3*t5*(rw*t2*t1172-rw*t7*t1170);
        double et62 = t1409;
        double et63 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et64 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et65 = t1410;
        double et66 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et67 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et68 = Itzz*(Iwyy*(-t1060+t1141+(t947*(t275-t415))/2.0+rw*t2*t8*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0))+rw*t3*t7*(t836-t1069+t190*t243*(t911+(t611*(t198-t363))/2.0)))+rw*t2*t1350-Iwyy*t8*(t961+t236*(t443+t673)+(t680*(t441+t527))/4.0+rw*t3*t7*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673))))-rw*t8*(Itzz*(t8*t321*(t815+t938+t1200-(t495*(t275-t415))/2.0)+t3*t10*t190*t1339)+Itzz*t3*t5*(rw*t2*t991+t3*t10*t190*t588+Itzz*Iwyy*t3*t8*t228*t243*t321));
        double et69 = Itzz*t10*(Itzz*t10*(t815+t938+t1200-(t495*(t275-t415))/2.0)+t3*t5*t87*t1168)+Itzz*t3*t5*(Iwyy*t1170+rw*t2*t1208);
        double et70 = t1409;
        double et71 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et72 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et73 = rw*t376*(Itzz*(Iwyy*t1362+Iwyy*t8*t1297-rw*t7*t1350)-t10*t87*(Itzz*t10*(-t817+t1198+(t516*(t275-t415))/2.0+t190*t272*(t347-t847))-t3*t5*t87*(t551+rw*t7*t942-rw*t2*t3*t802+(Iwyy*t190*t272*(t276-t491))/2.0))-rw*t8*(Itzz*(t8*t321*(-t817+t1198+(t516*(t275-t415))/2.0+t190*t272*(t347-t847))-t3*t10*t190*(-Iwyy*(t911+(t611*(t198-t363))/2.0)+Iwyy*t8*t757+rw*t7*t1247+rw*t7*t8*t870))-Itzz*t3*t5*(rw*t7*t991+t3*t10*t190*t596+Iwyy*t3*t8*t87*t228*t272*t321))+t3*t5*t87*(Iwyy*t1172+rw*t7*t1208));
        double et74 = 1.0/(et70+et71+et72);
        double et75 = (t997*(Itzz*t1406-t10*t87*t1342+rw*t72*t321*(t426+Itzz*(t347+t529-t772-t847-t950+t1049))+Itzz*t3*t5*(Itzz*t10*t1175-t3*t5*t87*t1173)))/(et56+et57+et58)+t626*t1421*(et59+et60+et61)+t1104*t1422*(t1232+t1396-rw*t3*t8*t10*t190*(t426+Itzz*(t347+t529-t772-t847-t950+t1049)))+(t835*(t10*t1396+rw*t8*(Itzz*t3*t5*t1301+t3*t74*t190*t1287)+Itzz*t3*t5*t1382))/(et62+et63+et64)+((Itzz*t1382+rw*t8*(Itzz*t1301+t30*t59*t74*t228*t370)+t10*t87*t1230)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et65+et66+et67)+rw*t368*t1421*(et68+et69);
        double et76 = et73*et74+rw*t3*t122*t1421*(t8*t321*(t1232+t1396)+t3*t10*t190*(Itzz*t1406-t10*t87*t1342+Itzz*t3*t5*(Itzz*t10*t1175-t3*t5*t87*t1173)));
        double et77 = t1410;
        double et78 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et79 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);

        double et80 = t1410;
        double et81 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et82 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et83 = t1409;
        double et84 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et85 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et86 = t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))-rw*t8*(t1384+t8*t321*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))));
        double et87 = rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))));
        double et88 = t1409;
        double et89 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et90 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et91 = t10*t87*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+Itzz*t3*t5*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))));
        double et92 = -Itzz*rw*t3*t5*t8*(t882+Iwyy*(t783-t862)+rw*t2*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))));
        double et93 = t1409;
        double et94 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et95 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et96 = t1409;
        double et97 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et98 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et99 = t1409;
        double et100 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et101 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et102 = ((t1386+Itzz*rw*t3*t5*t8*t1139)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et77+et78+et79);
        double et103 = (t1104*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))-Itzz*t8*t10*t59*t238*t1076))/(et80+et81+et82)+(t835*(et86+et87))/(et83+et84+et85);
        double et104 = (t626*(-Itzz*t10*t1381+t3*t5*t87*(t1001-t1222-t1322+rw*t2*t3*(t650+rw*t2*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+Itzz*rw*t3*t5*t8*(-t602+t1068+rw*t7*(t783-t862))))/(et88+et89+et90)+t997*t1422*(t10*t87*(-t896+t916-t1189+t1329+t236*(t347-t847)+rw*t2*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)))+Itzz*t3*t5*(-t725+t917+t988+t1165+t1188+t1285)+Itzz*rw*t3*t5*t72*t321*t1076)+(rw*t368*(et91+et92))/(et93+et94+et95);
        double et105 = (rw*t376*(t10*t87*(t1129-t1276-t1375+rw*t7*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+Itzz*t3*t5*t1394+Itzz*rw*t3*t5*t8*t1337))/(et96+et97+et98);
        double et106 = (rw*t3*t122*(Itzz*t10*(t1384+t8*t321*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))+Itzz*t3*t5*(t1364+t8*t321*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))))/(et99+et100+et101);

        d_dPH[index] = et19+et20+et21+et22+et23;
        
        index = index + num_threads;
        // index1 = (grid_size[0] > 1) ? index : 0;
        index2 = (grid_size[1] > 1) ? index : 0;
        index3 = (grid_size[2] > 1) ? index : 0;
        index4 = (grid_size[3] > 1) ? index : 0;
        index5 = (grid_size[4] > 1) ? index : 0;
        index6 = (grid_size[5] > 1) ? index : 0;
        index7 = (grid_size[6] > 1) ? index : 0;
        index8 = (grid_size[7] > 1) ? index : 0;

        uindex1 = (active_actions[0] == 1) ? index : 0;
        uindex2 = (active_actions[1] == 1) ? index : 0;
    }
}

__global__ void dyn4_mex_continuous(double* const d_dTH, // outputs
    double* const PH, double* const TH, double* const OM, double* const dPH, double* const dTH,
    double* const dPS, double* const dOM, double* const dNE, // input states
    double const* const in1, double const* const in2, // input actions
    const double mw, const double Iwxx, const double Iwyy, const double Iwzz, 
    const double mf, const double Ifxx, const double Ifyy, const double Ifzz, 
    const double mt, const double Itxx, const double Ityy, const double Itzz, 
    const double nw, const double nt, const double rw, const double rf, const double rt, const double fcoeff, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3]
                        * grid_size[4] * grid_size[5] * grid_size[6] * grid_size[7];
    double th, ps, om, ne, dph, dth, dps, dom, dne, Tneta, Tomega;
    ps = 0; ne = 0;

    // int index1 = (grid_size[0] > 1) ? index : 0;
    int index2 = (grid_size[1] > 1) ? index : 0;
    int index3 = (grid_size[2] > 1) ? index : 0;
    int index4 = (grid_size[3] > 1) ? index : 0;
    int index5 = (grid_size[4] > 1) ? index : 0;
    int index6 = (grid_size[5] > 1) ? index : 0;
    int index7 = (grid_size[6] > 1) ? index : 0;
    int index8 = (grid_size[7] > 1) ? index : 0;

    int uindex1 = (active_actions[0] == 1) ? index : 0;
    int uindex2 = (active_actions[1] == 1) ? index : 0;

    while (index < num_elements)
    {
        th = TH[index2];
        om = OM[index3];
        dph = dPH[index4];
        dth = dTH[index5];
        dps = dPS[index6];
        dom = dOM[index7];
        dne = dNE[index8];

        Tomega = nw*in1[uindex1];
        Tneta = nt*in2[uindex2];

        double t2 = cos(ph);
        double t3 = cos(th);
        double t4 = cos(ps);
        double t5 = cos(om);
        double t6 = cos(ne);
        double t7 = sin(ph);
        double t8 = sin(th);
        double t9 = sin(ps);
        double t10 = sin(om);
        double t11 = sin(ne);
        double t12 = mf+mt;
        double t13 = mf*rf;
        double t14 = mt*rt;
        double t15 = rf+rt;
        double t16 = Ifxx*2.0;
        double t17 = Ifxx*4.0;
        double t18 = Ifyy*2.0;
        double t19 = Ifyy*4.0;
        double t20 = Ifzz*2.0;
        double t21 = Ifzz*4.0;
        double t22 = Itxx*2.0;
        double t23 = Itxx*3.0;
        double t24 = Itxx*4.0;
        double t25 = Ityy*2.0;
        double t26 = Ityy*3.0;
        double t27 = Ityy*4.0;
        double t28 = Itzz*2.0;
        double t29 = Itzz*4.0;
        double t30 = Itzz*Itzz;
        double t31 = Iwxx*2.0;
        double t32 = Iwxx*4.0;
        double t33 = Iwyy*2.0;
        double t34 = Iwyy*4.0;
        double t35 = Iwzz*2.0;
        double t36 = Iwzz*4.0;
        double t37 = fcoeff*2.0;

        double t40 = mt*2.0;
        double t42 = rf*rf;

        double t44 = rw*rw;
        double t45 = Tomega*4.0;
        double t46 = th*2.0;
        double t47 = ps*2.0;
        double t48 = om*2.0;
        double t49 = ne*2.0;
        double t50 = dom*2.0;
        double t51 = dph*dph;
        double t52 = dth*dth;
        double t53 = dom*dom;
        double t79 = Ifyy*8.0;
        double t80 = -Ifzz;
        double t83 = -Ityy;
        double t87 = -Itzz;
        double t92 = Iwyy*8.0;
        double t93 = -Iwzz;
        double t109 = Ifxx/2.0;
        double t110 = Ifzz/2.0;
        double t111 = Itxx/4.0;
        double t112 = Ityy/4.0;
        double t113 = Itzz/2.0;
        double t114 = Iwxx/2.0;
        double t115 = Iwzz/2.0;
        double t54 = cos(t46);
        double t55 = cos(t47);
        double t56 = cos(t48);
        double t57 = cos(t49);
        double t58 = t2*t2;
        double t59 = t3*t3;
        double t60 = t4*t4;
        double t61 = t5*t5;
        double t62 = t5*t5*t5;
        double t63 = t6*t6;
        double t64 = t13*2.0;
        double t65 = t13*4.0;
        double t66 = t14*2.0;
        double t67 = sin(t46);
        double t68 = sin(t47);
        double t69 = sin(t48);
        double t70 = sin(t49);
        double t71 = t7*t7;
        double t72 = t8*t8;
        double t73 = t9*t9;
        double t74 = t10*t10;
        double t75 = t11*t11;
        double t76 = mw+t12;
        double t77 = -t16;
        double t78 = -t17;
        double t81 = -t20;
        double t82 = -t21;
        double t84 = -t25;
        double t85 = -t26;
        double t86 = -t27;
        double t88 = -t28;
        double t89 = -t29;
        double t90 = -t31;
        double t91 = -t32;
        double t94 = -t35;
        double t95 = -t36;
        double t96 = t8*dph;
        double t97 = -t45;
        double t98 = rf*t12;
        double t99 = mt*t15;
        double t100 = t2*t5;
        double t101 = t2*t10;
        double t102 = t5*t7;
        double t103 = t6*t8;
        double t104 = t2*dph*dps;
        double t105 = t7*t10;
        double t106 = t8*t11;
        double t107 = rf*t13;
        double t108 = t15*t15;
        double t116 = rf*t14*4.0;
        double t117 = rf*t14*6.0;
        double t118 = Ifyy*Iwyy*t8;
        double t120 = t15*t40;
        double t122 = -t52;
        double t124 = Itxx+t83;
        double t125 = Iwxx+t93;
        double t133 = rt*t14*3.0;
        double t134 = rt*t14*4.0;
        double t135 = rf*t14*8.0;
        double t140 = t3*t6*t10;
        double t141 = t51+t52;
        double t144 = t3*t10*t11;
        double t148 = t20+t28;
        double t155 = t12*2.0;
        double t156 = t12*3.0;
        double t157 = t2*t3*dph*dth*2.0;
        double t158 = t7*t8*t52;
        double t159 = t3*t7*dph*dth*2.0;
        double t179 = Iwyy*mt*t8*t42;
        double t180 = Iwyy*rt*t8*t14;
        double t181 = (Itxx*Iwyy*t8)/2.0;
        double t182 = (Ityy*Iwyy*t8)/2.0;
        double t183 = rf*t8*t14*t33;
        double t191 = Iwyy*rw*t3*t7*t8;
        double t219 = t12*t42*4.0;
        double t119 = t98*2.0;
        double t121 = t99*4.0;
        double t123 = t50*t96;
        double t126 = rw*t76;
        double t127 = t96+dps;
        double t128 = t96+dom;
        double t129 = t54*3.0;
        double t130 = rf*t64;
        double t131 = rf*t65;
        double t132 = rt*t66;

        double t137 = Itxx*t63;
        double t138 = Ityy*t75;
        double t139 = t8*t100;
        double t142 = t8*t101;
        double t143 = t8*t102;
        double t145 = t8*t105;
        double t146 = t15*t99;
        double t149 = t22*t63;
        double t150 = t27*t63;
        double t151 = t35*t60;
        double t152 = t24*t75;
        double t153 = t25*t75;
        double t154 = t31*t73;
        double t160 = Iwyy*t13*t101;
        double t161 = Iwyy*mt*rf*t101;
        double t162 = Iwyy*t14*t101;
        double t163 = Iwyy*t13*t105;
        double t164 = Iwyy*mt*rf*t105;
        double t165 = Iwyy*t14*t105;
        double t166 = t40*t108;
        double t168 = t51*t54;
        double t169 = t51*t59;
        double t171 = t107/2.0;
        double t172 = t22+t84;
        double t173 = t23+t85;
        double t174 = t24+t86;
        double t175 = t124*t124;
        double t176 = t31+t94;
        double t177 = t32+t95;
        double t178 = Iwyy*t8*t107;
        double t184 = t50+t96;
        double t186 = -t140;
        double t187 = -t157;
        double t189 = t14+t98;
        double t190 = t13+t99;
        double t203 = t56*t59*2.0;
        double t205 = t42*t155;
        double t206 = t42*t156;

        double t208 = Ityy*Iwyy*t11*t140;
        double t210 = t64+t120;
        double t212 = t57*t124;
        double t213 = t55*t125;
        double t214 = t68*t125;
        double t215 = t44*t58*t76;

        double t217 = t57*t182;
        double t218 = t44*t71*t76;
        double t223 = t2*t8*t141;
        double t225 = t61*t148;
        double t226 = t10*t57*t67*4.0;
        double t232 = Itxx*Iwyy*t8*t57*(-1.0/2.0);
        double t236 = t5*t6*t11*t124;
        double t244 = t103+t144;
        double t246 = t8*t70*t124;
        double t253 = t6*t11*t53*t124;
        double t254 = t4*t9*t52*t125;

        double t275 = t3*t10*t70*t124;
        double t287 = t5*t67*t70*t124;
        double t288 = t3*t69*t70*t124;

        double t333 = t6*t11*t61*t122*t124;
        double t336 = t3*t6*t11*t61*t87*t124;
        double t147 = -t129;
        double t167 = t15*t121;

        double t185 = -t139;
        double t188 = -t145;

        double t195 = Iwyy*t13*t143;
        double t196 = Iwyy*mt*rf*t143;
        double t197 = Iwyy*t14*t143;
        double t198 = t2*t8*t126;
        double t199 = t7*t8*t126;
        double t200 = t7*t127*dph;
        double t201 = t146/2.0;
        double t202 = -t168;
        double t204 = t190*t190;
        double t209 = t66+t119;
        double t211 = t65+t121;
        double t227 = t5*t189;
        double t228 = t5*t190;
        double t229 = t212*2.0;
        double t230 = t213*2.0;
        double t231 = t213*4.0;
        double t233 = t51+t53+t123;
        double t234 = t57*t173;
        double t235 = t55*t176;
        double t237 = t68*t176;
        double t239 = t107+t146;
        double t240 = Itxx*Iwyy*t11*t186;
        double t241 = Itzz+t212;
        double t242 = t101+t143;
        double t243 = t102+t142;

        double t247 = -t212;
        double t250 = t8*t218;
        double t256 = Itzz*Iwyy*t126*t142;

        double t258 = t246*4.0;
        double t259 = t8*t215;

        double t263 = t68*t177*dth*dps;
        double t267 = rw*t5*t210;
        double t274 = t106+t186;
        double t276 = t5*t246;

        double t279 = t3*t236*4.0;
        double t280 = g*t3*t10*t190*4.0;
        double t282 = rw*t8*t10*t190*4.0;
        double t283 = Iwyy*t87*t126*t145;
        double t286 = t3*t5*t70*t172;
        double t290 = t11*t124*t186;
        double t291 = t57*t61*t172;
        double t292 = t70*t124*t128;

        double t299 = t61*t70*t172*dth;
        double t301 = t288*2.0;
        double t305 = t4*t9*t125*t169;
        double t306 = t288*dph;
        double t307 = t275/2.0;
        double t311 = t16+t149+t153;
        double t312 = t10*t70*t96*t174*dth;
        double t315 = t10*t67*t70*t172;
        double t330 = t6*t11*t100*t124*t126;
        double t331 = t6*t11*t102*t124*t126;
        double t338 = t10*t54*t70*t174*dph;
        double t339 = t51*t287;
        double t346 = t61*t63*t75*t175;
        double t350 = rw*t30*t59*t62*t190;
        double t358 = t215*t236;
        double t368 = t104+t159+t223;
        double t370 = Iwyy+t215+t218;
        double t380 = Itxx+Ityy+t16+t81+t88+t130+t166;
        double t220 = Iwyy*t13*t185;
        double t221 = Iwyy*mt*rf*t185;
        double t222 = Iwyy*t14*t185;

        double t238 = rw*t228;
        double t248 = -t229;
        double t249 = -t231;
        double t251 = t234*2.0;
        double t252 = t235*2.0;
        double t255 = rw*t227*4.0;
        double t264 = -t234;
        double t266 = rw*t5*t209;
        double t268 = rw*t5*t211;
        double t269 = rw*t10*t209;
        double t270 = rw*t211*dth*dom;
        double t271 = Itzz+t247;
        double t272 = t100+t188;
        double t273 = t105+t185;

        double t281 = t235/4.0;
        double t284 = t28+t229;
        double t285 = t33+t230;
        double t289 = -t280;

        double t297 = t54*t239;
        double t298 = t56*t239;
        double t300 = t59*t237*dps;
        double t308 = t276/2.0;

        double t316 = rw*t10*t72*t211;
        double t317 = t3*t10*t241*4.0;
        double t318 = t126+t227;
        double t319 = -t307;
        double t321 = t126+t228;
        double t324 = t147+t203+1.0;
        double t327 = -t299;
        double t329 = rw*t8*t10*t51*t211;
        double t332 = t291/4.0;
        double t334 = -t305;
        double t337 = t53+t123+t202;
        double t343 = -t315;

        double t351 = -t338;
        double t353 = t74*t311;
        double t354 = t190*t242;
        double t356 = -t346;
        double t361 = Iwyy*t190*t243;
        double t366 = rw*t10*t190*t233;

        double t372 = rw*t5*t204*t243;
        double t374 = t228*t290;
        double t375 = rw*t2*t3*t190*t243;
        double t376 = t158+t187+t200;
        double t379 = rw*t7*t8*t190*t243;
        double t387 = t2*t126*t190*t243;
        double t398 = t56*t380;
        double t399 = Itxx+Ityy+t18+t130+t166+t247;
        double t401 = t279+t282;
        double t406 = Ifxx+t80+t87+t137+t138+t239;
        double t409 = t6*t11*t124*t228*t243;
        double t426 = t30*t59*t61*t370;
        double t428 = t51*t124*t244*t274;
        double t438 = t212+t380;
        double t458 = t17+t22+t25+t82+t89+t131+t167+t229;
        double t495 = t160+t161+t162+t195+t196+t197;
        double t663 = t118+t178+t179+t180+t181+t182+t183+t208+t217+t232+t240;
        double t260 = Iwyy+t238;
        double t265 = -t251;
        double t295 = t28+t248;
        double t296 = t34+t249;
        double t303 = t268/4.0;
        double t309 = t285*dps;
        double t314 = t5*t271*dth;

        double t325 = -t297;
        double t326 = -t298;
        double t328 = t67*t268*dph;
        double t335 = t3*t10*t271;
        double t340 = t5*t8*t284;
        double t345 = Iwyy*t8*t10*t87*t238;
        double t359 = g*t8*t318*4.0;
        double t360 = t70*t324;
        double t362 = Iwyy*t2*t3*t321;
        double t363 = t190*t273;
        double t364 = t5*t128*t284;
        double t365 = t3*t184*t268*dph;
        double t369 = t8*t361;
        double t371 = t214+t269;
        double t373 = -t366;
        double t378 = Iwyy+t213+t266;
        double t381 = rw*t59*t71*t321;
        double t382 = rw*t2*t3*t7*t8*t321;
        double t383 = Itzz*t59*t100*t321;
        double t384 = Itzz*t59*t102*t321;
        double t386 = rw*t5*t204*t272;
        double t388 = -t375;
        double t389 = rw*t2*t8*t190*t272;
        double t390 = rw*t3*t7*t190*t272;
        double t393 = t3*t58*t126*t321;
        double t394 = t3*t71*t126*t321;
        double t395 = t70*t124*t337;
        double t396 = t7*t126*t190*t272;
        double t404 = t199+t354;

        double t408 = t399*dom;
        double t410 = t19+t34+t150+t152+t255;
        double t411 = t398/4.0;
        double t412 = t3*t6*t11*t100*t124*t321;
        double t413 = rw*t3*t100*t190*t321;
        double t414 = t3*t6*t11*t102*t124*t321;
        double t415 = t8*t399;
        double t416 = rw*t3*t102*t190*t321;
        double t417 = t53*t401;
        double t419 = t19+t22+t25+t131+t167+t248;
        double t420 = t6*t11*t124*t228*t272;
        double t421 = t59*t101*t190*t321;

        double t423 = t59*t105*t190*t321;
        double t424 = (Iwyy*t399)/2.0;
        double t433 = t69*t406;

        double t453 = t3*t10*t438;
        double t454 = t8*t190*t243*t321;
        double t455 = t69*t438;
        double t461 = (t2*t126*t399)/2.0;
        double t462 = (t7*t126*t399)/2.0;
        double t468 = t3*t7*t190*t243*t321;
        double t470 = t8*t190*t272*t321;
        double t478 = t2*t3*t190*t272*t321;
        double t480 = t267+t399;
        double t482 = t69*t458*dth*dom;
        double t488 = t3*t56*t458*dph*dth;
        double t494 = t3*t56*t184*t438*dph;

        double t509 = (t2*t3*t321*t399)/2.0;
        double t511 = (t3*t7*t321*t399)/2.0;
        double t514 = t3*t56*t458*(t96-dom);
        double t516 = t163+t164+t165+t220+t221+t222;
        double t533 = t236*t495;
        double t563 = Itxx+Ityy+t19+t34+t77+t81+t88+t90+t94+t116+t132+t205+t235+t264;
        double t569 = t190*t272*t495;
        double t606 = (t399*t495)/2.0;
        double t680 = Itxx+Ityy+t16+t31+t35+t130+t148+t166+t235+t268+t291+t398;
        double t691 = t236*t663;
        double t717 = t190*t243*t663;
        double t720 = t2*t3*t321*t663;
        double t721 = t3*t7*t321*t663;
        double t723 = t190*t272*t663;
        double t302 = t8*t260;
        double t323 = t296*dps;
        double t341 = t335*dph;
        double t342 = -t314;
        double t344 = t5*t295*dom;
        double t347 = Iwyy*t72*t260;
        double t352 = t5*t8*t295*2.0;

        double t385 = t5*t11*t103*t124*t260;
        double t391 = t3*t371;

        double t397 = t3*t378*dph*dth;
        double t403 = t246+t335;
        double t405 = -t396;

        double t429 = -t414;
        double t430 = Itzz*t10*t404;
        double t431 = t419*dom;
        double t432 = -t424;

        double t435 = t288+t340;
        double t436 = t415/2.0;

        double t440 = t72*t410;
        double t445 = t3*t419*dph*dth;
        double t448 = t2*t126*t404;
        double t449 = -Itzz*t10*(t198-t363);
        double t450 = Iwyy*t10*t113*t415;
        double t456 = t236+t390;
        double t459 = t59*t433*2.0;

        double t465 = -t461;
        double t466 = -t462;
        double t471 = -t7*t126*(t198-t363);
        double t472 = -t124*dph*(t226-t360);
        double t473 = t236*t404;
        double t474 = t270+t395;
        double t475 = t238*t404;

        double t485 = t214+t433;
        double t489 = t236*(t198-t363);

        double t492 = (t8*t480)/2.0;
        double t493 = t309+t408;
        double t498 = t331+t413;
        double t499 = t3*t7*t321*t404;
        double t500 = t237+t455;

        double t506 = t190*t272*t404;

        double t513 = t455*(t52-t169);
        double t515 = -t190*t243*(t198-t363);

        double t518 = t393+t394;
        double t520 = -t5*(t246-t453);
        double t521 = t306+t327+t364;
        double t522 = Iwyy*t3*t5*t113*t321*t415;

        double t527 = (t7*t126*(t275-t415))/2.0;

        double t529 = rw*t2*t8*t516;
        double t540 = t151+t154+t225+t326+t353;
        double t543 = t236*(t275-t415)*(-1.0/2.0);
        double t546 = (t399*t404)/2.0;
        double t551 = t236*t516;

        double t573 = t3*t7*t321*(t275-t415)*(-1.0/2.0);

        double t577 = t421+t454;
        double t581 = t190*t243*t516;
        double t592 = t22+t25+t78+t79+t82+t89+t91+t92+t95+t134+t135+t219+t252+t265;
        double t600 = (t51*t67*t563)/2.0;

        double t619 = t468+t478;
        double t621 = t409+t509;
        double t622 = ((t198-t363)*(t275-t415))/2.0;
        double t634 = (t399*t516)/2.0;

        double t657 = t238*(t423-t470);

        double t671 = -rw*t3*t7*(t420-t511);
        double t676 = rw*t2*t8*(t420-t511);

        double t694 = ((t275-t415)*(t330-t416))/2.0;

        double t696 = (Itzz*t3*t5*t680)/4.0;
        double t703 = (t7*t126*t680)/4.0;
        double t712 = (Iwyy*t199*t680)/4.0;
        double t713 = (t215*t680)/4.0;

        double t725 = (Iwyy*t8*t238*t680)/4.0;
        double t741 = (t190*t243*t680)/4.0;
        double t744 = (t190*t272*t680)/4.0;
        double t759 = (t399*t680)/8.0;
        double t763 = (t404*t680)/4.0;

        double t792 = t663*(t275-t415)*(-1.0/2.0);

        double t810 = t680*(t275-t415)*(-1.0/8.0);
        double t825 = t109+t110+t111+t112+t113+t114+t115+t171+t201+t281+t303+t332+t381+t411;

        double t400 = t259+t302;
        double t402 = t391/2.0;
        double t418 = t403*dph*dom;
        double t441 = t190*t243*t302;
        double t442 = t435*dph;
        double t443 = t2*t3*t302*t321;
        double t444 = t3*t7*t302*t321;
        double t446 = t190*t272*t302;

        double t464 = -t459;
        double t469 = Iwyy*t456;

        double t486 = t10*t474*2.0;
        double t487 = t302*t404;
        double t491 = t3*t485;

        double t501 = t321*t436;
        double t502 = t3*t493*dph*2.0;

        double t510 = t292+t341+t342;
        double t512 = t8*t500;
        double t526 = t387+t405;
        double t530 = t521*dne*2.0;
        double t531 = Iwyy*t8*t518;
        double t532 = rw*t2*t518;
        double t534 = t383+t430;

        double t539 = rw*t3*t7*t518;
        double t548 = t384+t449;
        double t549 = t290+t492;
        double t552 = t388+t456;
        double t553 = t372+t466;
        double t555 = t236*t518;
        double t557 = t386+t465;
        double t558 = t59*t540*2.0;
        double t559 = t3*t10*t190*t518;
        double t575 = t190*t272*t498;

        double t598 = t448+t471;
        double t599 = t319+t389+t436;
        double t604 = t96*t592;
        double t611 = t391+t520;
        double t620 = t238*t577;
        double t626 = t254+t334+t373+t397+Tomega;
        double t631 = Iwyy*t621;
        double t642 = (t399*t518)/2.0;

        double t652 = rw*t2*t3*t619;
        double t654 = rw*t2*t8*t619;
        double t655 = rw*t3*t7*t619;
        double t656 = rw*t7*t8*t619;

        double t700 = t506+t515;
        double t714 = -Iwyy*(t499+t2*t3*t321*(t198-t363));
        double t749 = t302*(t499+t2*t3*t321*(t198-t363));
        double t802 = t569+t581;
        double t822 = t356+t759;
        double t824 = t412+t741;
        double t827 = Iwyy*t825;
        double t828 = t429+t744;
        double t832 = -t2*t126*(t346-t759);

        double t836 = ((t499+t2*t3*t321*(t198-t363))*(t275-t415))/2.0;
        double t840 = t7*t126*(t346-t759);
        double t842 = -rw*t2*t8*(t414-t744);
        double t867 = t302*(t346-t759);
        double t898 = t551+t721;
        double t930 = t634+t723;

        double t467 = t236*t400;
        double t476 = t250+t400;
        double t497 = t491/2.0;
        double t519 = t510*dne*4.0;

        double t535 = -t530;
        double t536 = rw*t2*t526;
        double t538 = t8*t532;
        double t541 = Iwyy*t534;
        double t542 = dth*(t344-t442)*(-1.0/2.0);
        double t545 = rw*t3*t7*t526;

        double t550 = Iwyy*t548;
        double t556 = Itzz*t10*t549;
        double t560 = rw*t7*t553;
        double t561 = t286+t512;

        double t565 = rw*t2*(t276-t491)*(-1.0/2.0);
        double t566 = t2*t126*t549;
        double t567 = t7*t126*t549;
        double t568 = rw*t2*t557;

        double t571 = t238*t534;
        double t574 = t287+t316+t464;
        double t578 = t3*t5*t87*t552;

        double t583 = t238*t548;
        double t584 = t236*t549;
        double t589 = t8*t321*t526;
        double t590 = t374+t501;
        double t593 = t236*(t276-t491)*(-1.0/2.0);
        double t601 = t302*t549;
        double t603 = Iwyy*t599;
        double t607 = t362+t532;
        double t608 = rw*t2*t598;
        double t616 = rw*t3*t7*t598;
        double t618 = t190*t243*t549;
        double t623 = t2*t3*t321*t549;
        double t624 = t3*t7*t321*t549;
        double t627 = t190*t272*t549;
        double t630 = -t620;
        double t633 = t190*t243*(t276-t491)*(-1.0/2.0);

        double t636 = t3*t7*t321*(t276-t491)*(-1.0/2.0);
        double t637 = t361*(t276-t491)*(-1.0/2.0);
        double t643 = (t2*t3*t321*(t276-t491))/2.0;

        double t649 = t190*t272*(t276-t491)*(-1.0/2.0);
        double t650 = t8*t631;

        double t658 = t236*t598;
        double t662 = t3*t5*t113*t611;
        double t666 = (t2*t126*t611)/2.0;
        double t667 = (t7*t126*t611)/2.0;
        double t668 = (t399*t534)/2.0;
        double t672 = Itzz*t3*t74*t190*t598;
        double t678 = t399*(t276-t491)*(-1.0/4.0);
        double t679 = (t399*t548)/2.0;
        double t681 = (t236*t611)/2.0;
        double t682 = (t238*t611)/2.0;
        double t685 = (t3*t10*t190*t611)/2.0;
        double t689 = -rw*t2*(t446+(t2*t126*(t275-t415))/2.0);

        double t693 = t379+t599;
        double t698 = (t190*t243*t611)/2.0;
        double t708 = (t190*t272*t611)/2.0;

        double t716 = rw*t2*t700;
        double t719 = (t399*t598)/2.0;
        double t728 = ((t275-t415)*(t276-t491))/4.0;

        double t731 = (t399*t611)/4.0;
        double t733 = (t404*t611)/2.0;
        double t735 = -Iwyy*(t382-t402+(t5*(t246-t453))/2.0);

        double t747 = t526*(t276-t491)*(-1.0/2.0);

        double t756 = (t400*t680)/4.0;
        double t760 = (t598*(t275-t415))/2.0;
        double t767 = t323+t431+t604;
        double t813 = (t549*t611)/2.0;

        double t815 = rw*t2*t8*t802;
        double t816 = rw*t3*t7*t802;
        double t817 = rw*t7*t8*t802;
        double t820 = (t611*(t276-t491))/4.0;
        double t830 = (t526*t680)/4.0;
        double t834 = t117+t133+t206+t325+t343+t440+t558;

        double t843 = ((t441+t527)*(t276-t491))/2.0;
        double t849 = t238*t824;
        double t850 = (t549*t680)/4.0;
        double t852 = ((t276-t491)*(t446+(t2*t126*(t275-t415))/2.0))/2.0;
        double t864 = t539+t703;
        double t870 = t531+t714;
        double t896 = (t611*t663)/2.0;
        double t905 = Itzz*t3*t5*(t696-t10*t113*(t276-t491));
        double t916 = rw*t2*t8*t898;
        double t917 = rw*t3*t7*t898;

        double t929 = t680*(t446+(t2*t126*(t275-t415))/2.0)*(-1.0/4.0);
        double t935 = rw*t2*t8*t930;
        double t936 = rw*t3*t7*t930;
        double t971 = -Iwyy*t8*(-t575+t642+t190*t243*(t330-t416));

        double t976 = -Itzz*t3*t5*(-t575+t642+t190*t243*(t330-t416));

        double t979 = -rw*t3*t7*(-t575+t642+t190*t243*(t330-t416));
        double t980 = rw*t2*(-t575+t642+t190*t243*(t330-t416));
        double t982 = t655+t824;

        double t989 = t652+t828;
        double t990 = t671+t822;

        double t1086 = -rw*t2*(t832+t238*(t414-t744));
        double t1090 = Iwyy*t8*(t832+t238*(t414-t744));
        double t1092 = rw*t2*(t832+t238*(t414-t744));

        double t544 = t8*t536;
        double t579 = t52*t561;
        double t582 = t574*dom;

        double t588 = t283+t541;

        double t594 = Itzz*t10*t399*t476*(-1.0/2.0);
        double t595 = t191+t565;
        double t596 = t256+t550;
        double t597 = Iwyy*t590;
        double t609 = t8*t603;
        double t610 = t361+t536;

        double t613 = t2*t126*t590;
        double t614 = t7*t126*t590;
        double t615 = t8*t608;
        double t628 = -t618;
        double t629 = t336+t556;

        double t639 = rw*t2*t3*t607;
        double t640 = rw*t7*t8*t607;
        double t648 = rw*t7*(t331-t545);

        double t665 = t302*t590;
        double t673 = -t667;

        double t701 = Itzz*t10*t693;

        double t710 = -Iwyy*(t308-t497+rw*t3*t7*(t198-t363));
        double t722 = t8*t716;

        double t736 = -rw*t7*(t475-t567);
        double t742 = rw*t2*(t566-t238*(t198-t363));
        double t743 = rw*t7*t8*t87*(t475-t567);
        double t748 = t385+t682;
        double t751 = t473+t623;
        double t754 = t473+t633;
        double t757 = t444+t666;
        double t758 = t3*t74*t87*t190*(t475-t567);
        double t761 = t495+t608;

        double t765 = Itzz*t3*t74*t190*(t566-t238*(t198-t363));

        double t774 = t3*t767;
        double t775 = t489+t649;
        double t778 = t190*t243*(t566-t238*(t198-t363));

        double t797 = t236*(t489-t624);
        double t799 = t3*t10*t190*(t489-t624);
        double t800 = t559+t589;

        double t835 = t253+t333+t418+t428+t542+Tneta;
        double t841 = rw*t7*(t616-(t7*t126*(t276-t491))/2.0);
        double t846 = -t843;
        double t847 = (Iwyy*t834)/4.0;
        double t853 = (Itzz*t10*t834)/4.0;
        double t854 = ((t275-t415)*(t475-t567))/2.0;

        double t859 = (t2*t126*t834)/4.0;
        double t860 = (t7*t126*t834)/4.0;

        double t865 = ((t566-t238*(t198-t363))*(t275-t415))/2.0;
        double t873 = t8*t321*(t627-(t399*(t198-t363))/2.0);
        double t875 = rw*t7*t864;
        double t876 = t584+t678;
        double t879 = (t236*t834)/4.0;
        double t880 = (t238*t834)/4.0;
        double t887 = t432+t560+t568;
        double t888 = rw*t2*t8*t870;
        double t889 = rw*t3*t7*t870;

        double t903 = t238*(t685-(t8*t321*(t275-t415))/2.0);
        double t904 = t573+t708;
        double t906 = ((t275-t415)*(t489-t624))/2.0;
        double t907 = (t190*t243*t834)/4.0;
        double t910 = (t2*t3*t321*t834)/4.0;
        double t911 = (t3*t7*t321*t834)/4.0;
        double t913 = (t190*t272*t834)/4.0;

        double t942 = t658+t747;
        double t947 = t643+t763;
        double t952 = t555+t830;
        double t959 = -t7*t126*(t636+(t680*(t198-t363))/4.0);
        double t985 = t593+t850;
        double t993 = t681+t810;

        double t995 = rw*t7*t8*t982;
        double t996 = Iwyy*t990;
        double t997 = t97+t289+t312+t339+t445+t488+t513+t519;
        double t1000 = t3*t5*t87*t982;
        double t1004 = Itzz*t3*t5*t989;
        double t1011 = t399*(t636+(t680*(t198-t363))/4.0)*(-1.0/2.0);
        double t1015 = (t680*t834)/1.6e+1;
        double t1027 = -Itzz*t10*(-t654+t698+(t2*t3*t321*(t275-t415))/2.0);
        double t1030 = -rw*t2*t3*(-t654+t698+(t2*t3*t321*(t275-t415))/2.0);
        double t1031 = t543+t676+t731;
        double t1033 = ((t636+(t680*(t198-t363))/4.0)*(t275-t415))/2.0;
        double t1073 = t840+t849;
        double t1074 = t631+t980;
        double t1091 = t8*t1086;
        double t1104 = t263+t359+t365+t482+t486+t494+t502+t535+t600;
        double t587 = -t582;
        double t602 = t8*t597;

        double t625 = t2*t126*t595;
        double t638 = Iwyy*t629;
        double t644 = rw*t2*t3*t610;
        double t646 = rw*t7*t8*t610;

        double t661 = t10*t87*t629;
        double t664 = t7*t126*t629;
        double t669 = t2*t126*t629;
        double t697 = t190*t243*t629;
        double t706 = t190*t272*t629;
        double t711 = (t399*t588)/2.0;
        double t718 = (t399*t596)/2.0;
        double t745 = Itzz*t8*t742;
        double t766 = rw*t2*t757;
        double t768 = t2*t126*t751;

        double t771 = rw*t3*t7*t757;
        double t772 = rw*t7*t8*t761;
        double t776 = rw*t2*t3*t761;
        double t777 = Iwyy*t8*(t443+t673);
        double t781 = t236*t751;

        double t783 = t3*t10*t190*t751;
        double t784 = t238*t754;
        double t785 = rw*t2*t8*t775;
        double t788 = t236*t757;
        double t793 = t302*t751;
        double t798 = t190*t243*t748;

        double t804 = t190*t272*t748;
        double t806 = t238*t775;
        double t807 = t546+t628;
        double t809 = t190*t272*t751;
        double t812 = t190*t243*t757;
        double t818 = Itzz*t3*t5*t800;

        double t829 = t441+t527+t544;
        double t845 = (t399*t757)/2.0;

        double t866 = -t860;
        double t871 = t614+t630;
        double t881 = -t880;
        double t884 = t613+t657;

        double t899 = t2*t126*t876;
        double t900 = t7*t126*t876;
        double t902 = t578+t701;
        double t908 = -t903;

        double t912 = t3*t5*t10*t30*t887;
        double t914 = -t910;
        double t920 = rw*t3*t7*t904;
        double t922 = t302*t876;
        double t923 = t443+t538+t673;
        double t938 = t190*t243*(t347-t847);
        double t945 = rw*t2*t942;

        double t949 = Iwyy*t947;
        double t950 = rw*t2*(t859-t302*(t198-t363));
        double t953 = t2*t126*t947;
        double t955 = t328+t351+t514+t774;
        double t958 = t236*(t859-t302*(t198-t363));

        double t961 = rw*t2*t8*t952;
        double t962 = rw*t7*t8*t952;
        double t963 = t238*t947;
        double t965 = t3*t5*t87*t952;
        double t970 = t302*t947;
        double t975 = t190*t243*(t859-t302*(t198-t363));
        double t984 = t190*t272*t947;
        double t987 = Itzz*t985;
        double t988 = Iwyy*t985;
        double t1001 = t8*t996;
        double t1002 = t2*t126*t985;
        double t1003 = t7*t126*t985;
        double t1007 = (t399*(t859-t302*(t198-t363)))/2.0;
        double t1008 = (t399*t947)/2.0;
        double t1013 = t662+t853;
        double t1014 = t238*t993;

        double t1018 = t190*t243*t985;
        double t1019 = t190*t272*t985;
        double t1023 = t656+t904;
        double t1024 = t622+t913;
        double t1036 = Iwyy*t1031;
        double t1042 = t728+t879;
        double t1043 = t985*(t275-t415)*(-1.0/2.0);
        double t1058 = -t7*t126*(t911+(t611*(t198-t363))/2.0);
        double t1064 = -t236*(t911+(t611*(t198-t363))/2.0);
        double t1066 = t799+t873;
        double t1071 = t813+t879;
        double t1076 = t663+t736+t742;
        double t1078 = Iwyy*t8*t1073;
        double t1083 = rw*t2*t3*t1074;
        double t1084 = rw*t7*t8*t1074;
        double t1095 = t399*(t911+(t611*(t198-t363))/2.0)*(-1.0/2.0);
        double t1101 = -rw*t7*(-t722+t907+(t404*(t275-t415))/2.0);
        double t1106 = -rw*t2*(-t719+t778+t190*t272*(t475-t567));
        double t1120 = t820+t1015;
        double t1170 = t1000+t1027;
        double t1173 = t639+t713+t827+t875;
        double t1182 = t979+t1073;

        double t1204 = t10*t87*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567)));
        double t1206 = -rw*t7*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567)));
        double t659 = -t644;
        double t674 = -t664;
        double t675 = t37+t300+t587;
        double t702 = t345+t638;
        double t707 = -t697;
        double t838 = rw*t7*t829;
        double t848 = -t845;
        double t862 = t8*t321*t807;
        double t863 = t583+t669;
        double t882 = Iwyy*t8*t871;
        double t883 = rw*t7*t871;
        double t892 = Itzz*t3*t5*t871;
        double t894 = Iwyy*t8*t884;
        double t895 = rw*t2*t884;
        double t901 = Itzz*t3*t5*t884;

        double t924 = t238*t902;
        double t925 = rw*t7*t923;
        double t933 = t679+t706;
        double t943 = t487+t866;
        double t948 = t8*t945;
        double t957 = t955*dth;
        double t991 = t672+t818;

        double t1006 = -t1003;
        double t1009 = -t1008;
        double t1010 = t601+t881;
        double t1016 = Itzz*t10*t1013;
        double t1020 = t756+t771;
        double t1028 = rw*t2*t1024;
        double t1029 = Itzz*t10*t1023;
        double t1034 = rw*t3*t7*t1024;
        double t1035 = t665+t908;
        double t1040 = t8*t1036;
        double t1044 = t238*t1042;
        double t1046 = t733+t914;
        double t1067 = Iwyy*t1066;
        double t1068 = rw*t2*t1066;

        double t1080 = t2*t126*t1071;
        double t1081 = t7*t126*t1071;
        double t1082 = -rw*t2*(t806-t899);
        double t1093 = t236*t1071;
        double t1094 = t661+t987;
        double t1096 = t190*t243*t1071;
        double t1097 = t190*t272*t1071;
        double t1107 = t8*t1106;

        double t1110 = t852+t958;
        double t1111 = -rw*t2*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624));
        double t1113 = -rw*t3*t7*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624));
        double t1117 = -t3*t10*t190*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624));
        double t1125 = t87*t1120;
        double t1126 = -Itzz*t10*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673));
        double t1127 = -Iwyy*t8*(t798+(t399*(t443+t673))/2.0+(t498*(t275-t415))/2.0);
        double t1130 = -rw*t2*t3*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673));
        double t1134 = t2*t126*t1120;
        double t1135 = t7*t126*t1120;
        double t1141 = t190*t243*t1120;
        double t1142 = t190*t272*t1120;
        double t1146 = ((t275-t415)*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624)))/2.0;
        double t1148 = (t399*t1120)/2.0;
        double t1152 = -rw*t3*t7*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624));
        double t1156 = t953+t959;
        double t1159 = ((t276-t491)*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673)))/2.0;

        double t1163 = Itzz*(t1002-t238*(t636+(t680*(t198-t363))/4.0));
        double t1165 = rw*t2*(t1002-t238*(t636+(t680*(t198-t363))/4.0));
        double t1166 = t625+t710+t776+t841;
        double t1168 = t533+t637+t816+t945;
        double t1180 = rw*t7*t8*(t650+rw*t2*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)));
        double t1183 = rw*t7*t1182;
        double t1196 = -rw*t7*t8*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0));
        double t1211 = -rw*t7*(t798+(t399*(t443+t673))/2.0+t8*t980+(t498*(t275-t415))/2.0);
        double t1213 = -rw*t2*t3*(t798+(t399*(t443+t673))/2.0+t8*t980+(t498*(t275-t415))/2.0);
        double t1219 = t797+t1011+t1019;
        double t1233 = -rw*t2*t3*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)));
        double t1244 = -rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)));
        double t1291 = -Iwyy*(t1120+rw*t2*t8*(t636+(t680*(t198-t363))/4.0)+rw*t3*t7*(t911+(t611*(t198-t363))/2.0));
        double t1300 = -rw*t7*(t961+t236*(t443+t673)+(t680*(t441+t527))/4.0+rw*t3*t7*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673)));
        double t1323 = t842+t920+t993+t995+t1030;

        double t686 = t675*dph*2.0;
        double t737 = t3*t10*t190*t702;
        double t750 = t190*t243*t702;
        double t752 = t190*t272*t702;
        double t844 = -t838;
        double t855 = t571+t674;
        double t872 = rw*t2*t863;
        double t886 = t10*t87*t863;
        double t928 = -t925;
        double t931 = t668+t707;
        double t940 = rw*t2*t933;
        double t954 = t236*t943;
        double t973 = t190*t272*t943;
        double t981 = t358+t469+t648+t659;
        double t999 = (t399*t943)/2.0;
        double t1012 = Itzz*t1010;
        double t1022 = Iwyy*t8*t1020;
        double t1025 = t594+t924;
        double t1032 = t190*t243*t1010;
        double t1038 = t190*t272*t1010;
        double t1045 = t615+t943;
        double t1047 = -t1044;
        double t1048 = Iwyy*t1046;
        double t1052 = t758+t892;
        double t1053 = t2*t126*t1046;

        double t1057 = t765+t901;
        double t1060 = t236*t1046;
        double t1061 = t238*t1046;
        double t1069 = t190*t272*t1046;
        double t1070 = -t1068;
        double t1087 = -t1081;
        double t1088 = (t399*t1046)/2.0;
        double t1089 = t8*t1082;
        double t1099 = t302*t1094;
        double t1112 = t8*t1111;
        double t1114 = rw*t2*t1110;
        double t1116 = Itzz*t1113;
        double t1123 = t694+t804+t848;
        double t1138 = -t1135;
        double t1139 = t597+t883+t895;
        double t1149 = -t1148;

        double t1157 = rw*t2*t1156;
        double t1160 = t963+t1006;
        double t1167 = Itzz*t3*t5*t1166;
        double t1171 = rw*t7*t8*t1168;
        double t1172 = t1004+t1029;
        double t1184 = t533+t720+t1111;
        double t1185 = t8*t1183;
        double t1186 = -t1183;
        double t1208 = t965+t1126;
        double t1210 = t781+t1009+t1018;
        double t1216 = rw*t2*(t1134-t302*(t636+(t680*(t198-t363))/4.0));
        double t1221 = Iwyy*t1219;
        double t1222 = rw*t2*t1219;
        double t1225 = t7*t126*t1219;
        double t1227 = t302*t1219;
        double t1240 = t785+t1034+t1042;
        double t1246 = t8*t1244;
        double t1269 = t190*t243*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624));

        double t1271 = -rw*t7*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567)));
        double t1272 = t905+t1016+t1125;
        double t1275 = t906+t1095+t1097;
        double t1294 = -rw*t3*t7*(Itzz*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624))+t3*t5*t87*t863);
        double t1297 = t788+t929+t962+t1130;
        double t1303 = t1033+t1064+t1142;
        double t1315 = -t7*t126*(-t1093+t1148+(t985*(t275-t415))/2.0);
        double t1316 = t7*t126*(-t1093+t1148+(t985*(t275-t415))/2.0);
        double t1324 = Iwyy*t1323;

        double t868 = Itzz*t10*t855;
        double t869 = rw*t7*t855;
        double t877 = Itzz*t3*t5*t855;
        double t932 = rw*t7*t931;
        double t941 = -t940;
        double t998 = t3*t5*t87*t981;
        double t1039 = -t1032;
        double t1041 = t3*t5*t87*t1025;
        double t1049 = rw*t7*t1045;
        double t1055 = rw*t7*t1052;
        double t1063 = rw*t2*t1057;
        double t1115 = Itzz*t1112;
        double t1128 = rw*t2*t1123;
        double t1129 = Iwyy*t8*t1123;
        double t1137 = rw*t3*t7*t1123;
        double t1140 = t603+t646+t689+t844;
        double t1158 = t8*t1157;

        double t1164 = t87*t1160;
        double t1169 = t350+t743+t745+t1012;
        double t1175 = t640+t735+t766+t928;
        double t1188 = rw*t2*t3*t1184;
        double t1189 = rw*t7*t8*t1184;
        double t1194 = t760+t973+t975;
        double t1212 = Iwyy*t1210;
        double t1215 = t886+t1163;
        double t1217 = t2*t126*t1210;
        double t1218 = -t1216;
        double t1223 = t8*t1222;
        double t1224 = t302*t1210;
        double t1237 = t865+t1007+t1038;
        double t1241 = Iwyy*t1240;
        double t1247 = t749+t1053+t1058;
        double t1253 = t793+t1061+t1087;
        double t1268 = rw*t2*t3*(t836-t1069+t190*t243*(t911+(t611*(t198-t363))/2.0));

        double t1276 = Iwyy*t1275;
        double t1277 = t238*t1272;
        double t1278 = rw*t2*t1275;
        double t1279 = rw*t3*t7*t1275;
        double t1282 = t718+t752+t1206;
        double t1285 = -rw*t7*(t1160+rw*t3*t7*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624)));
        double t1296 = rw*t7*t8*(t712+t889-t949-t1157);
        double t1299 = rw*t2*t1297;
        double t1304 = -t238*(-t1060+t1141+(t947*(t275-t415))/2.0);
        double t1305 = t238*t1303;
        double t1311 = t1043+t1093+t1149;
        double t1321 = t1152+t1210;
        double t1336 = -rw*t7*(-t1088+t1096+(t751*(t275-t415))/2.0+rw*t2*t8*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)));
        double t1337 = t894+t1067+t1271;
        double t1344 = t996+t1083+t1092+t1186;
        double t874 = -t869;
        double t1133 = -t1128;
        double t1144 = Itzz*t10*t1140;
        double t1154 = t8*t321*t1140;
        double t1155 = t450+t932+t941;
        double t1181 = t3*t10*t190*t1175;

        double t1197 = rw*t2*t1194;
        double t1198 = rw*t7*t1194;
        double t1220 = rw*t2*t8*t1215;
        double t1236 = t854+t999+t1039;
        double t1238 = rw*t2*t1237;
        double t1239 = rw*t3*t7*t1237;
        double t1242 = (t680*t1194)/4.0;
        double t1248 = rw*t2*t1247;
        double t1249 = rw*t3*t7*t1247;
        double t1250 = t236*t1247;
        double t1252 = (t399*t1247)/2.0;
        double t1254 = Itzz*t1253;
        double t1258 = t522+t737+t1055+t1063;
        double t1262 = t190*t272*t1253;
        double t1280 = -t1278;
        double t1287 = -Itzz*(t347+t529-t772-t847-t950+t1049);
        double t1306 = t868+t1116+t1164;
        double t1322 = rw*t7*t1321;
        double t1325 = t867+t1014+t1091+t1137;
        double t1327 = t1112+t1253;
        double t1345 = Itzz*t3*t5*t1344;
        double t1352 = -Itzz*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)));
        double t1356 = rw*t7*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)));
        double t1362 = t1196+t1268+t1303;
        double t1364 = t3*t10*t190*(-t725+t917+t988+t1165+t1188+t1285);
        double t1368 = t1223+t1279+t1311;
        double t1382 = t1299+t1300+t1324;
        double t1153 = t702+t872+t874;
        double t1199 = t3*t1198;
        double t1200 = -t1197;

        double t1230 = t998+t1144;
        double t1251 = -t1250;
        double t1260 = t3*t5*t87*t1258;
        double t1263 = -t1262;
        double t1288 = t10*t1287;
        double t1301 = t1154+t1181;
        double t1307 = rw*t7*t8*t1306;
        double t1317 = t1107+t1236;
        double t1326 = Iwyy*t8*t1325;
        double t1329 = rw*t7*t1327;
        double t1338 = t877+t1115+t1254;
        double t1339 = t777+t888+t1048+t1248;
        double t1346 = t1036+t1084+t1133+t1211;
        double t1348 = t970+t1138+t1158+t1249;
        double t1360 = t1185+t1213+t1325;
        double t1367 = t1204+t1352;
        double t1369 = Iwyy*t1368;
        double t1381 = t1040+t1180+t1280+t1336;
        double t1394 = t1090+t1221+t1233+t1356;
        double t1232 = Itzz*t3*t5*t1230;
        double t1313 = rw*t2*t3*(t815+t938+t1200-(t495*(t275-t415))/2.0);
        double t1318 = rw*t7*t1317;
        double t1319 = rw*t2*t3*t1317;

        double t1331 = t846+t948+t954+t1199;
        double t1340 = rw*t2*t3*t1338;
        double t1341 = rw*t2*t3*t1339;
        double t1342 = t1167+t1288;
        double t1347 = -t10*t87*t1346;
        double t1349 = rw*t7*t1348;
        double t1350 = t1159+t1242+t1251;
        double t1372 = t1146+t1252+t1263+t1269;
        double t1384 = -t3*t10*t190*(-t896+t916-t1189+t1329+t236*(t347-t847)+rw*t2*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)));

        double t1332 = rw*t7*t1331;
        double t1374 = rw*t2*t1372;
        double t1375 = rw*t7*t1372;
        double t1386 = t1345+t1347;
        double t1389 = t922+t1047+t1089+t1239+t1246+t1319;
        double t1405 = -Itzz*(t1384+t8*t321*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))));
        double t1406 = t1022+t1218+t1291+t1296+t1341+t1349;
        double t1333 = -t1332;
        double t1376 = t3*t1375;
        double t1390 = Itzz*t1389;
        double t1408 = t1260+t1405;
        double t1396 = Itzz*(t1114+t1171+t1241+t1313+t1333-Iwyy*t8*(t467+rw*t3*t7*(t446+(t2*t126*(t275-t415))/2.0)));
        double t1409 = rw*t8*t1408;
        double t1414 = -rw*t7*(-t1224+t1315+t1376+t238*(-t1060+t1141+(t947*(t275-t415))/2.0)+rw*t2*t8*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0))));
        double t1410 = t1409*4.0;

        double et1 = t1409;
        double et2 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et3 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double t1421 = -1.0/(et1+et2+et3);
        double et4 = t1410;
        double et5 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et6 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double t1422 = -1.0/(et4+et5+et6);

        double et7 = t1410;
        double et8 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et9 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et10 = t1409;
        double et11 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et12 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et13 = t1409;
        double et14 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et15 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et16 = t1409;
        double et17 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et18 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et19 = t1422*(Itzz*t1344+t30*t74*t887+Itzz*rw*t8*t1139)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352)))+t835*t1421*(t1345-Itzz*t10*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+rw*t8*t1258);
        double et20 = t1104*t1422*(t912+Itzz*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))-Itzz*rw*t3*t8*t10*t190*t1076)+(t997*(Itzz*(-t725+t917+t988+t1165+t1188+t1285)+t10*t87*t1153+Itzz*rw*t72*t321*t1076))/(et7+et8+et9);
        double et21 = (t626*(Itzz*(t1001-t1222-t1322+rw*t2*t3*(t650+rw*t2*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))-Itzz*t10*t1155+rw*t8*t87*(-t602+t1068+rw*t7*(t783-t862))))/(et10+et11+et12);
        double et22 = (rw*t368*(-Itzz*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+t10*t87*(-t711+t750+rw*t2*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567))))+Itzz*rw*t8*(t882+Iwyy*(t783-t862)+rw*t2*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))))))/(et13+et14+et15)+rw*t376*t1421*(Itzz*t1394+t10*t87*t1282+Itzz*rw*t8*t1337);
        double et23 = (rw*t3*t52*(Itzz*(t1364+t8*t321*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t10*t1258))/(et16+et17+et18);
        double et24 = t1410;
        double et25 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et26 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et27 = t1409;
        double et28 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et29 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et30 = t1409;
        double et31 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et32 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et33 = t1409;
        double et34 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et35 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et36 = t1409;
        double et37 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et38 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et39 = ((t912-Itzz*t1346)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et24+et25+et26)+(t835*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+t3*t5*t87*t1346))/(et27+et28+et29)+t1104*t1422*(Itzz*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))-t30*t59*t61*t887)+(t626*(Itzz*t1381+Itzz*t3*t5*t1155))/(et30+et31+et32);
        double et40 = t997*t1422*(Itzz*(-t896+t916-t1189+t1329+t236*(t347-t847)+rw*t2*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)))+t3*t5*t87*t1153)+(rw*t368*(Itzz*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))-t3*t5*t87*(-t711+t750+rw*t2*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567))))))/(et33+et34+et35);
        double et41 = (rw*t376*(Itzz*(t1129-t1276-t1375+rw*t7*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+t3*t5*t87*t1282))/(et36+et37+et38)+rw*t3*t52*t1408*t1421;
        double et42 = t1410;
        double et43 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et44 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et45 = t1409;
        double et46 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et47 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et48 = t1409;
        double et49 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et50 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et51 = t87*(-t1093+t1148+(t985*(t275-t415))/2.0)+Itzz*t10*((t399*t1013)/2.0+(t629*(t275-t415))/2.0)-rw*t8*(Itzz*((t590*t834)/4.0-t549*(t685-(t8*t321*(t275-t415))/2.0))+Itzz*t8*t1070-Itzz*t3*t5*(t3*t10*t190*t629+t3*t5*t113*t321*t415)+rw*t7*t8*t87*(t783-t862))-rw*t2*t3*(Itzz*(-t1088+t1096+(t751*(t275-t415))/2.0)-Itzz*t3*t5*t931+Itzz*rw*t2*t8*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))+rw*t3*t7*(Itzz*t1275-t3*t5*t87*t933);
        double et52 = Itzz*t3*t5*(t236*t629+(t399*(t696-t10*t113*(t276-t491)))/2.0)+rw*t7*t8*(Itzz*t1152+Itzz*t1210+Itzz*t10*t931)+rw*t2*t8*(Itzz*t1219+t10*t87*t933);
        double et53 = t997*t1422*(t1099+t1220+t1277+t1294+t1307+t1340+rw*t72*t321*t1169)+((t87*t1360+Itzz*t10*t1025+rw*t8*(Itzz*t1035+Itzz*t8*t883+Itzz*t8*t895-rw*t30*t59*t61*t74*t204))*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et42+et43+et44)+t1104*t1422*(t1041+t1390+rw*t3*t8*t10*t190*t1169)+(t835*(t10*t1390+rw*t8*(t8*t1055+t8*t1063+Itzz*t3*t5*t1035+t3*t74*t190*t1012)+t3*t5*t87*t1360))/(et45+et46+et47)+(t626*(et51+et52))/(et48+et49+et50);
        double et54 = rw*t376*t1421*(t87*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+Itzz*t10*(Itzz*t10*t1237+t3*t5*t87*t1123)+rw*t8*(t87*(t8*t321*t1237+t3*t10*t190*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)))+Itzz*t3*t5*t1057+Itzz*rw*t7*t8*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))))+rw*t7*t8*t1367+rw*t2*t3*(Itzz*t1372-Itzz*t3*t5*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567))))-Itzz*t3*t5*(Itzz*t10*(t806-t899)+t3*t5*t87*(t832+t238*(t414-t744))));
        double et55 = rw*t368*t1421*(Itzz*(t1224+t1304+t1316)-t10*t87*(Itzz*t10*t1236+Itzz*t3*t5*(t798+(t399*(t443+t673))/2.0+(t498*(t275-t415))/2.0))+rw*t8*(-Itzz*(t8*t321*t1236-t3*t10*t190*t1253)+t3*t5*t87*t1052+Itzz*rw*t2*t8*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))))-t3*t5*t87*(Itzz*t10*(t784-t900)-t3*t5*t87*t1073)+rw*t2*t8*t1367-rw*t3*t7*(Itzz*t1372-Itzz*t3*t5*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567)))))+rw*t3*t122*t1421*(t8*t321*(t1041+t1390)-t3*t10*t190*(t1099+t1220+t1277+t1294+t1307+t1340));
        double et56 = t1410;
        double et57 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et58 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et59 = -Itzz*(-t8*t1324+rw*t2*t1362+rw*t7*(-t1060+t1141+(t947*(t275-t415))/2.0+rw*t2*t8*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0))+rw*t3*t7*(t836-t1069+t190*t243*(t911+(t611*(t198-t363))/2.0))));
        double et60 = rw*t8*(Itzz*(t8*t321*(t609-t1028+t1101+rw*t7*t8*(t369-t716))-t3*t10*t190*(Iwyy*t8*(t382-t402+(t5*(t246-t453))/2.0)-rw*t2*(t911+(t611*(t198-t363))/2.0)+rw*t7*(-t733+t910+rw*t2*t8*(t499+t2*t3*t321*(t198-t363)))-rw*t7*t8*(t8*t362+rw*t2*(t499+t2*t3*t321*(t198-t363)))))-t3*t5*t87*(rw*t7*(t3*t10*t190*t534+Itzz*t3*t8*t228*t243*t321)-rw*t2*(t3*t10*t190*t548+t3*t8*t87*t228*t272*t321)+Itzz*Iwyy*t3*t8*t74*t190));
        double et61 = t10*t87*(Itzz*t10*(t609-t1028+t1101+rw*t7*t8*(t369-t716))-Itzz*t3*t5*(t8*t469+rw*t7*(t754+rw*t3*t7*t700)+rw*t2*t775-rw*t2*t3*(t369-t716)))+Itzz*t3*t5*(rw*t2*t1172-rw*t7*t1170);
        double et62 = t1409;
        double et63 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et64 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et65 = t1410;
        double et66 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et67 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et68 = Itzz*(Iwyy*(-t1060+t1141+(t947*(t275-t415))/2.0+rw*t2*t8*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0))+rw*t3*t7*(t836-t1069+t190*t243*(t911+(t611*(t198-t363))/2.0)))+rw*t2*t1350-Iwyy*t8*(t961+t236*(t443+t673)+(t680*(t441+t527))/4.0+rw*t3*t7*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673))))-rw*t8*(Itzz*(t8*t321*(t815+t938+t1200-(t495*(t275-t415))/2.0)+t3*t10*t190*t1339)+Itzz*t3*t5*(rw*t2*t991+t3*t10*t190*t588+Itzz*Iwyy*t3*t8*t228*t243*t321));
        double et69 = Itzz*t10*(Itzz*t10*(t815+t938+t1200-(t495*(t275-t415))/2.0)+t3*t5*t87*t1168)+Itzz*t3*t5*(Iwyy*t1170+rw*t2*t1208);
        double et70 = t1409;
        double et71 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et72 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et73 = rw*t376*(Itzz*(Iwyy*t1362+Iwyy*t8*t1297-rw*t7*t1350)-t10*t87*(Itzz*t10*(-t817+t1198+(t516*(t275-t415))/2.0+t190*t272*(t347-t847))-t3*t5*t87*(t551+rw*t7*t942-rw*t2*t3*t802+(Iwyy*t190*t272*(t276-t491))/2.0))-rw*t8*(Itzz*(t8*t321*(-t817+t1198+(t516*(t275-t415))/2.0+t190*t272*(t347-t847))-t3*t10*t190*(-Iwyy*(t911+(t611*(t198-t363))/2.0)+Iwyy*t8*t757+rw*t7*t1247+rw*t7*t8*t870))-Itzz*t3*t5*(rw*t7*t991+t3*t10*t190*t596+Iwyy*t3*t8*t87*t228*t272*t321))+t3*t5*t87*(Iwyy*t1172+rw*t7*t1208));
        double et74 = 1.0/(et70+et71+et72);
        double et75 = (t997*(Itzz*t1406-t10*t87*t1342+rw*t72*t321*(t426+Itzz*(t347+t529-t772-t847-t950+t1049))+Itzz*t3*t5*(Itzz*t10*t1175-t3*t5*t87*t1173)))/(et56+et57+et58)+t626*t1421*(et59+et60+et61)+t1104*t1422*(t1232+t1396-rw*t3*t8*t10*t190*(t426+Itzz*(t347+t529-t772-t847-t950+t1049)))+(t835*(t10*t1396+rw*t8*(Itzz*t3*t5*t1301+t3*t74*t190*t1287)+Itzz*t3*t5*t1382))/(et62+et63+et64)+((Itzz*t1382+rw*t8*(Itzz*t1301+t30*t59*t74*t228*t370)+t10*t87*t1230)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et65+et66+et67)+rw*t368*t1421*(et68+et69);
        double et76 = et73*et74+rw*t3*t122*t1421*(t8*t321*(t1232+t1396)+t3*t10*t190*(Itzz*t1406-t10*t87*t1342+Itzz*t3*t5*(Itzz*t10*t1175-t3*t5*t87*t1173)));
        double et77 = t1410;
        double et78 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et79 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);

        double et80 = t1410;
        double et81 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et82 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et83 = t1409;
        double et84 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et85 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et86 = t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))-rw*t8*(t1384+t8*t321*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))));
        double et87 = rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))));
        double et88 = t1409;
        double et89 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et90 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et91 = t10*t87*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+Itzz*t3*t5*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))));
        double et92 = -Itzz*rw*t3*t5*t8*(t882+Iwyy*(t783-t862)+rw*t2*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))));
        double et93 = t1409;
        double et94 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et95 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et96 = t1409;
        double et97 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et98 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et99 = t1409;
        double et100 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et101 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et102 = ((t1386+Itzz*rw*t3*t5*t8*t1139)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et77+et78+et79);
        double et103 = (t1104*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))-Itzz*t8*t10*t59*t238*t1076))/(et80+et81+et82)+(t835*(et86+et87))/(et83+et84+et85);
        double et104 = (t626*(-Itzz*t10*t1381+t3*t5*t87*(t1001-t1222-t1322+rw*t2*t3*(t650+rw*t2*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+Itzz*rw*t3*t5*t8*(-t602+t1068+rw*t7*(t783-t862))))/(et88+et89+et90)+t997*t1422*(t10*t87*(-t896+t916-t1189+t1329+t236*(t347-t847)+rw*t2*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)))+Itzz*t3*t5*(-t725+t917+t988+t1165+t1188+t1285)+Itzz*rw*t3*t5*t72*t321*t1076)+(rw*t368*(et91+et92))/(et93+et94+et95);
        double et105 = (rw*t376*(t10*t87*(t1129-t1276-t1375+rw*t7*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+Itzz*t3*t5*t1394+Itzz*rw*t3*t5*t8*t1337))/(et96+et97+et98);
        double et106 = (rw*t3*t122*(Itzz*t10*(t1384+t8*t321*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))+Itzz*t3*t5*(t1364+t8*t321*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))))/(et99+et100+et101);

        d_dTH[index] = et39+et40+et41;
        
        index = index + num_threads;
        // index1 = (grid_size[0] > 1) ? index : 0;
        index2 = (grid_size[1] > 1) ? index : 0;
        index3 = (grid_size[2] > 1) ? index : 0;
        index4 = (grid_size[3] > 1) ? index : 0;
        index5 = (grid_size[4] > 1) ? index : 0;
        index6 = (grid_size[5] > 1) ? index : 0;
        index7 = (grid_size[6] > 1) ? index : 0;
        index8 = (grid_size[7] > 1) ? index : 0;

        uindex1 = (active_actions[0] == 1) ? index : 0;
        uindex2 = (active_actions[1] == 1) ? index : 0;
    }
}

__global__ void dyn5_mex_continuous(double* const d_dPS, // outputs
    double* const PH, double* const TH, double* const OM, double* const dPH, double* const dTH,
    double* const dPS, double* const dOM, double* const dNE, // input states
    double const* const in1, double const* const in2, // input actions
    const double mw, const double Iwxx, const double Iwyy, const double Iwzz, 
    const double mf, const double Ifxx, const double Ifyy, const double Ifzz, 
    const double mt, const double Itxx, const double Ityy, const double Itzz, 
    const double nw, const double nt, const double rw, const double rf, const double rt, const double fcoeff, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3]
                        * grid_size[4] * grid_size[5] * grid_size[6] * grid_size[7];
    double th, ps, om, ne, dph, dth, dps, dom, dne, Tneta, Tomega;
    ps = 0; ne = 0;

    // int index1 = (grid_size[0] > 1) ? index : 0;
    int index2 = (grid_size[1] > 1) ? index : 0;
    int index3 = (grid_size[2] > 1) ? index : 0;
    int index4 = (grid_size[3] > 1) ? index : 0;
    int index5 = (grid_size[4] > 1) ? index : 0;
    int index6 = (grid_size[5] > 1) ? index : 0;
    int index7 = (grid_size[6] > 1) ? index : 0;
    int index8 = (grid_size[7] > 1) ? index : 0;

    int uindex1 = (active_actions[0] == 1) ? index : 0;
    int uindex2 = (active_actions[1] == 1) ? index : 0;

    while (index < num_elements)
    {
        th = TH[index2];
        om = OM[index3];
        dph = dPH[index4];
        dth = dTH[index5];
        dps = dPS[index6];
        dom = dOM[index7];
        dne = dNE[index8];

        Tomega = nw*in1[uindex1];
        Tneta = nt*in2[uindex2];

        double t2 = cos(ph);
        double t3 = cos(th);
        double t4 = cos(ps);
        double t5 = cos(om);
        double t6 = cos(ne);
        double t7 = sin(ph);
        double t8 = sin(th);
        double t9 = sin(ps);
        double t10 = sin(om);
        double t11 = sin(ne);
        double t12 = mf+mt;
        double t13 = mf*rf;
        double t14 = mt*rt;
        double t15 = rf+rt;
        double t16 = Ifxx*2.0;
        double t17 = Ifxx*4.0;
        double t18 = Ifyy*2.0;
        double t19 = Ifyy*4.0;
        double t20 = Ifzz*2.0;
        double t21 = Ifzz*4.0;
        double t22 = Itxx*2.0;
        double t23 = Itxx*3.0;
        double t24 = Itxx*4.0;
        double t25 = Ityy*2.0;
        double t26 = Ityy*3.0;
        double t27 = Ityy*4.0;
        double t28 = Itzz*2.0;
        double t29 = Itzz*4.0;
        double t30 = Itzz*Itzz;
        double t31 = Iwxx*2.0;
        double t32 = Iwxx*4.0;
        double t33 = Iwyy*2.0;
        double t34 = Iwyy*4.0;
        double t35 = Iwzz*2.0;
        double t36 = Iwzz*4.0;
        double t37 = fcoeff*2.0;

        double t40 = mt*2.0;
        double t42 = rf*rf;

        double t44 = rw*rw;
        double t45 = Tomega*4.0;
        double t46 = th*2.0;
        double t47 = ps*2.0;
        double t48 = om*2.0;
        double t49 = ne*2.0;
        double t50 = dom*2.0;
        double t51 = dph*dph;
        double t52 = dth*dth;
        double t53 = dom*dom;
        double t79 = Ifyy*8.0;
        double t80 = -Ifzz;
        double t83 = -Ityy;
        double t87 = -Itzz;
        double t92 = Iwyy*8.0;
        double t93 = -Iwzz;
        double t109 = Ifxx/2.0;
        double t110 = Ifzz/2.0;
        double t111 = Itxx/4.0;
        double t112 = Ityy/4.0;
        double t113 = Itzz/2.0;
        double t114 = Iwxx/2.0;
        double t115 = Iwzz/2.0;
        double t54 = cos(t46);
        double t55 = cos(t47);
        double t56 = cos(t48);
        double t57 = cos(t49);
        double t58 = t2*t2;
        double t59 = t3*t3;
        double t60 = t4*t4;
        double t61 = t5*t5;
        double t62 = t5*t5*t5;
        double t63 = t6*t6;
        double t64 = t13*2.0;
        double t65 = t13*4.0;
        double t66 = t14*2.0;
        double t67 = sin(t46);
        double t68 = sin(t47);
        double t69 = sin(t48);
        double t70 = sin(t49);
        double t71 = t7*t7;
        double t72 = t8*t8;
        double t73 = t9*t9;
        double t74 = t10*t10;
        double t75 = t11*t11;
        double t76 = mw+t12;
        double t77 = -t16;
        double t78 = -t17;
        double t81 = -t20;
        double t82 = -t21;
        double t84 = -t25;
        double t85 = -t26;
        double t86 = -t27;
        double t88 = -t28;
        double t89 = -t29;
        double t90 = -t31;
        double t91 = -t32;
        double t94 = -t35;
        double t95 = -t36;
        double t96 = t8*dph;
        double t97 = -t45;
        double t98 = rf*t12;
        double t99 = mt*t15;
        double t100 = t2*t5;
        double t101 = t2*t10;
        double t102 = t5*t7;
        double t103 = t6*t8;
        double t104 = t2*dph*dps;
        double t105 = t7*t10;
        double t106 = t8*t11;
        double t107 = rf*t13;
        double t108 = t15*t15;
        double t116 = rf*t14*4.0;
        double t117 = rf*t14*6.0;
        double t118 = Ifyy*Iwyy*t8;
        double t120 = t15*t40;
        double t122 = -t52;
        double t124 = Itxx+t83;
        double t125 = Iwxx+t93;
        double t133 = rt*t14*3.0;
        double t134 = rt*t14*4.0;
        double t135 = rf*t14*8.0;
        double t140 = t3*t6*t10;
        double t141 = t51+t52;
        double t144 = t3*t10*t11;
        double t148 = t20+t28;
        double t155 = t12*2.0;
        double t156 = t12*3.0;
        double t157 = t2*t3*dph*dth*2.0;
        double t158 = t7*t8*t52;
        double t159 = t3*t7*dph*dth*2.0;
        double t179 = Iwyy*mt*t8*t42;
        double t180 = Iwyy*rt*t8*t14;
        double t181 = (Itxx*Iwyy*t8)/2.0;
        double t182 = (Ityy*Iwyy*t8)/2.0;
        double t183 = rf*t8*t14*t33;
        double t191 = Iwyy*rw*t3*t7*t8;
        double t219 = t12*t42*4.0;
        double t119 = t98*2.0;
        double t121 = t99*4.0;
        double t123 = t50*t96;
        double t126 = rw*t76;
        double t127 = t96+dps;
        double t128 = t96+dom;
        double t129 = t54*3.0;
        double t130 = rf*t64;
        double t131 = rf*t65;
        double t132 = rt*t66;

        double t137 = Itxx*t63;
        double t138 = Ityy*t75;
        double t139 = t8*t100;
        double t142 = t8*t101;
        double t143 = t8*t102;
        double t145 = t8*t105;
        double t146 = t15*t99;
        double t149 = t22*t63;
        double t150 = t27*t63;
        double t151 = t35*t60;
        double t152 = t24*t75;
        double t153 = t25*t75;
        double t154 = t31*t73;
        double t160 = Iwyy*t13*t101;
        double t161 = Iwyy*mt*rf*t101;
        double t162 = Iwyy*t14*t101;
        double t163 = Iwyy*t13*t105;
        double t164 = Iwyy*mt*rf*t105;
        double t165 = Iwyy*t14*t105;
        double t166 = t40*t108;
        double t168 = t51*t54;
        double t169 = t51*t59;
        double t171 = t107/2.0;
        double t172 = t22+t84;
        double t173 = t23+t85;
        double t174 = t24+t86;
        double t175 = t124*t124;
        double t176 = t31+t94;
        double t177 = t32+t95;
        double t178 = Iwyy*t8*t107;
        double t184 = t50+t96;
        double t186 = -t140;
        double t187 = -t157;
        double t189 = t14+t98;
        double t190 = t13+t99;
        double t203 = t56*t59*2.0;
        double t205 = t42*t155;
        double t206 = t42*t156;

        double t208 = Ityy*Iwyy*t11*t140;
        double t210 = t64+t120;
        double t212 = t57*t124;
        double t213 = t55*t125;
        double t214 = t68*t125;
        double t215 = t44*t58*t76;

        double t217 = t57*t182;
        double t218 = t44*t71*t76;
        double t223 = t2*t8*t141;
        double t225 = t61*t148;
        double t226 = t10*t57*t67*4.0;
        double t232 = Itxx*Iwyy*t8*t57*(-1.0/2.0);
        double t236 = t5*t6*t11*t124;
        double t244 = t103+t144;
        double t246 = t8*t70*t124;
        double t253 = t6*t11*t53*t124;
        double t254 = t4*t9*t52*t125;

        double t275 = t3*t10*t70*t124;
        double t287 = t5*t67*t70*t124;
        double t288 = t3*t69*t70*t124;

        double t333 = t6*t11*t61*t122*t124;
        double t336 = t3*t6*t11*t61*t87*t124;
        double t147 = -t129;
        double t167 = t15*t121;

        double t185 = -t139;
        double t188 = -t145;

        double t195 = Iwyy*t13*t143;
        double t196 = Iwyy*mt*rf*t143;
        double t197 = Iwyy*t14*t143;
        double t198 = t2*t8*t126;
        double t199 = t7*t8*t126;
        double t200 = t7*t127*dph;
        double t201 = t146/2.0;
        double t202 = -t168;
        double t204 = t190*t190;
        double t209 = t66+t119;
        double t211 = t65+t121;
        double t227 = t5*t189;
        double t228 = t5*t190;
        double t229 = t212*2.0;
        double t230 = t213*2.0;
        double t231 = t213*4.0;
        double t233 = t51+t53+t123;
        double t234 = t57*t173;
        double t235 = t55*t176;
        double t237 = t68*t176;
        double t239 = t107+t146;
        double t240 = Itxx*Iwyy*t11*t186;
        double t241 = Itzz+t212;
        double t242 = t101+t143;
        double t243 = t102+t142;

        double t247 = -t212;
        double t250 = t8*t218;
        double t256 = Itzz*Iwyy*t126*t142;

        double t258 = t246*4.0;
        double t259 = t8*t215;

        double t263 = t68*t177*dth*dps;
        double t267 = rw*t5*t210;
        double t274 = t106+t186;
        double t276 = t5*t246;

        double t279 = t3*t236*4.0;
        double t280 = g*t3*t10*t190*4.0;
        double t282 = rw*t8*t10*t190*4.0;
        double t283 = Iwyy*t87*t126*t145;
        double t286 = t3*t5*t70*t172;
        double t290 = t11*t124*t186;
        double t291 = t57*t61*t172;
        double t292 = t70*t124*t128;

        double t299 = t61*t70*t172*dth;
        double t301 = t288*2.0;
        double t305 = t4*t9*t125*t169;
        double t306 = t288*dph;
        double t307 = t275/2.0;
        double t311 = t16+t149+t153;
        double t312 = t10*t70*t96*t174*dth;
        double t315 = t10*t67*t70*t172;
        double t330 = t6*t11*t100*t124*t126;
        double t331 = t6*t11*t102*t124*t126;
        double t338 = t10*t54*t70*t174*dph;
        double t339 = t51*t287;
        double t346 = t61*t63*t75*t175;
        double t350 = rw*t30*t59*t62*t190;
        double t358 = t215*t236;
        double t368 = t104+t159+t223;
        double t370 = Iwyy+t215+t218;
        double t380 = Itxx+Ityy+t16+t81+t88+t130+t166;
        double t220 = Iwyy*t13*t185;
        double t221 = Iwyy*mt*rf*t185;
        double t222 = Iwyy*t14*t185;

        double t238 = rw*t228;
        double t248 = -t229;
        double t249 = -t231;
        double t251 = t234*2.0;
        double t252 = t235*2.0;
        double t255 = rw*t227*4.0;
        double t264 = -t234;
        double t266 = rw*t5*t209;
        double t268 = rw*t5*t211;
        double t269 = rw*t10*t209;
        double t270 = rw*t211*dth*dom;
        double t271 = Itzz+t247;
        double t272 = t100+t188;
        double t273 = t105+t185;

        double t281 = t235/4.0;
        double t284 = t28+t229;
        double t285 = t33+t230;
        double t289 = -t280;

        double t297 = t54*t239;
        double t298 = t56*t239;
        double t300 = t59*t237*dps;
        double t308 = t276/2.0;

        double t316 = rw*t10*t72*t211;
        double t317 = t3*t10*t241*4.0;
        double t318 = t126+t227;
        double t319 = -t307;
        double t321 = t126+t228;
        double t324 = t147+t203+1.0;
        double t327 = -t299;
        double t329 = rw*t8*t10*t51*t211;
        double t332 = t291/4.0;
        double t334 = -t305;
        double t337 = t53+t123+t202;
        double t343 = -t315;

        double t351 = -t338;
        double t353 = t74*t311;
        double t354 = t190*t242;
        double t356 = -t346;
        double t361 = Iwyy*t190*t243;
        double t366 = rw*t10*t190*t233;

        double t372 = rw*t5*t204*t243;
        double t374 = t228*t290;
        double t375 = rw*t2*t3*t190*t243;
        double t376 = t158+t187+t200;
        double t379 = rw*t7*t8*t190*t243;
        double t387 = t2*t126*t190*t243;
        double t398 = t56*t380;
        double t399 = Itxx+Ityy+t18+t130+t166+t247;
        double t401 = t279+t282;
        double t406 = Ifxx+t80+t87+t137+t138+t239;
        double t409 = t6*t11*t124*t228*t243;
        double t426 = t30*t59*t61*t370;
        double t428 = t51*t124*t244*t274;
        double t438 = t212+t380;
        double t458 = t17+t22+t25+t82+t89+t131+t167+t229;
        double t495 = t160+t161+t162+t195+t196+t197;
        double t663 = t118+t178+t179+t180+t181+t182+t183+t208+t217+t232+t240;
        double t260 = Iwyy+t238;
        double t265 = -t251;
        double t295 = t28+t248;
        double t296 = t34+t249;
        double t303 = t268/4.0;
        double t309 = t285*dps;
        double t314 = t5*t271*dth;

        double t325 = -t297;
        double t326 = -t298;
        double t328 = t67*t268*dph;
        double t335 = t3*t10*t271;
        double t340 = t5*t8*t284;
        double t345 = Iwyy*t8*t10*t87*t238;
        double t359 = g*t8*t318*4.0;
        double t360 = t70*t324;
        double t362 = Iwyy*t2*t3*t321;
        double t363 = t190*t273;
        double t364 = t5*t128*t284;
        double t365 = t3*t184*t268*dph;
        double t369 = t8*t361;
        double t371 = t214+t269;
        double t373 = -t366;
        double t378 = Iwyy+t213+t266;
        double t381 = rw*t59*t71*t321;
        double t382 = rw*t2*t3*t7*t8*t321;
        double t383 = Itzz*t59*t100*t321;
        double t384 = Itzz*t59*t102*t321;
        double t386 = rw*t5*t204*t272;
        double t388 = -t375;
        double t389 = rw*t2*t8*t190*t272;
        double t390 = rw*t3*t7*t190*t272;
        double t393 = t3*t58*t126*t321;
        double t394 = t3*t71*t126*t321;
        double t395 = t70*t124*t337;
        double t396 = t7*t126*t190*t272;
        double t404 = t199+t354;

        double t408 = t399*dom;
        double t410 = t19+t34+t150+t152+t255;
        double t411 = t398/4.0;
        double t412 = t3*t6*t11*t100*t124*t321;
        double t413 = rw*t3*t100*t190*t321;
        double t414 = t3*t6*t11*t102*t124*t321;
        double t415 = t8*t399;
        double t416 = rw*t3*t102*t190*t321;
        double t417 = t53*t401;
        double t419 = t19+t22+t25+t131+t167+t248;
        double t420 = t6*t11*t124*t228*t272;
        double t421 = t59*t101*t190*t321;

        double t423 = t59*t105*t190*t321;
        double t424 = (Iwyy*t399)/2.0;
        double t433 = t69*t406;

        double t453 = t3*t10*t438;
        double t454 = t8*t190*t243*t321;
        double t455 = t69*t438;
        double t461 = (t2*t126*t399)/2.0;
        double t462 = (t7*t126*t399)/2.0;
        double t468 = t3*t7*t190*t243*t321;
        double t470 = t8*t190*t272*t321;
        double t478 = t2*t3*t190*t272*t321;
        double t480 = t267+t399;
        double t482 = t69*t458*dth*dom;
        double t488 = t3*t56*t458*dph*dth;
        double t494 = t3*t56*t184*t438*dph;

        double t509 = (t2*t3*t321*t399)/2.0;
        double t511 = (t3*t7*t321*t399)/2.0;
        double t514 = t3*t56*t458*(t96-dom);
        double t516 = t163+t164+t165+t220+t221+t222;
        double t533 = t236*t495;
        double t563 = Itxx+Ityy+t19+t34+t77+t81+t88+t90+t94+t116+t132+t205+t235+t264;
        double t569 = t190*t272*t495;
        double t606 = (t399*t495)/2.0;
        double t680 = Itxx+Ityy+t16+t31+t35+t130+t148+t166+t235+t268+t291+t398;
        double t691 = t236*t663;
        double t717 = t190*t243*t663;
        double t720 = t2*t3*t321*t663;
        double t721 = t3*t7*t321*t663;
        double t723 = t190*t272*t663;
        double t302 = t8*t260;
        double t323 = t296*dps;
        double t341 = t335*dph;
        double t342 = -t314;
        double t344 = t5*t295*dom;
        double t347 = Iwyy*t72*t260;
        double t352 = t5*t8*t295*2.0;

        double t385 = t5*t11*t103*t124*t260;
        double t391 = t3*t371;

        double t397 = t3*t378*dph*dth;
        double t403 = t246+t335;
        double t405 = -t396;

        double t429 = -t414;
        double t430 = Itzz*t10*t404;
        double t431 = t419*dom;
        double t432 = -t424;

        double t435 = t288+t340;
        double t436 = t415/2.0;

        double t440 = t72*t410;
        double t445 = t3*t419*dph*dth;
        double t448 = t2*t126*t404;
        double t449 = -Itzz*t10*(t198-t363);
        double t450 = Iwyy*t10*t113*t415;
        double t456 = t236+t390;
        double t459 = t59*t433*2.0;

        double t465 = -t461;
        double t466 = -t462;
        double t471 = -t7*t126*(t198-t363);
        double t472 = -t124*dph*(t226-t360);
        double t473 = t236*t404;
        double t474 = t270+t395;
        double t475 = t238*t404;

        double t485 = t214+t433;
        double t489 = t236*(t198-t363);

        double t492 = (t8*t480)/2.0;
        double t493 = t309+t408;
        double t498 = t331+t413;
        double t499 = t3*t7*t321*t404;
        double t500 = t237+t455;

        double t506 = t190*t272*t404;

        double t513 = t455*(t52-t169);
        double t515 = -t190*t243*(t198-t363);

        double t518 = t393+t394;
        double t520 = -t5*(t246-t453);
        double t521 = t306+t327+t364;
        double t522 = Iwyy*t3*t5*t113*t321*t415;

        double t527 = (t7*t126*(t275-t415))/2.0;

        double t529 = rw*t2*t8*t516;
        double t540 = t151+t154+t225+t326+t353;
        double t543 = t236*(t275-t415)*(-1.0/2.0);
        double t546 = (t399*t404)/2.0;
        double t551 = t236*t516;

        double t573 = t3*t7*t321*(t275-t415)*(-1.0/2.0);

        double t577 = t421+t454;
        double t581 = t190*t243*t516;
        double t592 = t22+t25+t78+t79+t82+t89+t91+t92+t95+t134+t135+t219+t252+t265;
        double t600 = (t51*t67*t563)/2.0;

        double t619 = t468+t478;
        double t621 = t409+t509;
        double t622 = ((t198-t363)*(t275-t415))/2.0;
        double t634 = (t399*t516)/2.0;

        double t657 = t238*(t423-t470);

        double t671 = -rw*t3*t7*(t420-t511);
        double t676 = rw*t2*t8*(t420-t511);

        double t694 = ((t275-t415)*(t330-t416))/2.0;

        double t696 = (Itzz*t3*t5*t680)/4.0;
        double t703 = (t7*t126*t680)/4.0;
        double t712 = (Iwyy*t199*t680)/4.0;
        double t713 = (t215*t680)/4.0;

        double t725 = (Iwyy*t8*t238*t680)/4.0;
        double t741 = (t190*t243*t680)/4.0;
        double t744 = (t190*t272*t680)/4.0;
        double t759 = (t399*t680)/8.0;
        double t763 = (t404*t680)/4.0;

        double t792 = t663*(t275-t415)*(-1.0/2.0);

        double t810 = t680*(t275-t415)*(-1.0/8.0);
        double t825 = t109+t110+t111+t112+t113+t114+t115+t171+t201+t281+t303+t332+t381+t411;

        double t400 = t259+t302;
        double t402 = t391/2.0;
        double t418 = t403*dph*dom;
        double t441 = t190*t243*t302;
        double t442 = t435*dph;
        double t443 = t2*t3*t302*t321;
        double t444 = t3*t7*t302*t321;
        double t446 = t190*t272*t302;

        double t464 = -t459;
        double t469 = Iwyy*t456;

        double t486 = t10*t474*2.0;
        double t487 = t302*t404;
        double t491 = t3*t485;

        double t501 = t321*t436;
        double t502 = t3*t493*dph*2.0;

        double t510 = t292+t341+t342;
        double t512 = t8*t500;
        double t526 = t387+t405;
        double t530 = t521*dne*2.0;
        double t531 = Iwyy*t8*t518;
        double t532 = rw*t2*t518;
        double t534 = t383+t430;

        double t539 = rw*t3*t7*t518;
        double t548 = t384+t449;
        double t549 = t290+t492;
        double t552 = t388+t456;
        double t553 = t372+t466;
        double t555 = t236*t518;
        double t557 = t386+t465;
        double t558 = t59*t540*2.0;
        double t559 = t3*t10*t190*t518;
        double t575 = t190*t272*t498;

        double t598 = t448+t471;
        double t599 = t319+t389+t436;
        double t604 = t96*t592;
        double t611 = t391+t520;
        double t620 = t238*t577;
        double t626 = t254+t334+t373+t397+Tomega;
        double t631 = Iwyy*t621;
        double t642 = (t399*t518)/2.0;

        double t652 = rw*t2*t3*t619;
        double t654 = rw*t2*t8*t619;
        double t655 = rw*t3*t7*t619;
        double t656 = rw*t7*t8*t619;

        double t700 = t506+t515;
        double t714 = -Iwyy*(t499+t2*t3*t321*(t198-t363));
        double t749 = t302*(t499+t2*t3*t321*(t198-t363));
        double t802 = t569+t581;
        double t822 = t356+t759;
        double t824 = t412+t741;
        double t827 = Iwyy*t825;
        double t828 = t429+t744;
        double t832 = -t2*t126*(t346-t759);

        double t836 = ((t499+t2*t3*t321*(t198-t363))*(t275-t415))/2.0;
        double t840 = t7*t126*(t346-t759);
        double t842 = -rw*t2*t8*(t414-t744);
        double t867 = t302*(t346-t759);
        double t898 = t551+t721;
        double t930 = t634+t723;

        double t467 = t236*t400;
        double t476 = t250+t400;
        double t497 = t491/2.0;
        double t519 = t510*dne*4.0;

        double t535 = -t530;
        double t536 = rw*t2*t526;
        double t538 = t8*t532;
        double t541 = Iwyy*t534;
        double t542 = dth*(t344-t442)*(-1.0/2.0);
        double t545 = rw*t3*t7*t526;

        double t550 = Iwyy*t548;
        double t556 = Itzz*t10*t549;
        double t560 = rw*t7*t553;
        double t561 = t286+t512;

        double t565 = rw*t2*(t276-t491)*(-1.0/2.0);
        double t566 = t2*t126*t549;
        double t567 = t7*t126*t549;
        double t568 = rw*t2*t557;

        double t571 = t238*t534;
        double t574 = t287+t316+t464;
        double t578 = t3*t5*t87*t552;

        double t583 = t238*t548;
        double t584 = t236*t549;
        double t589 = t8*t321*t526;
        double t590 = t374+t501;
        double t593 = t236*(t276-t491)*(-1.0/2.0);
        double t601 = t302*t549;
        double t603 = Iwyy*t599;
        double t607 = t362+t532;
        double t608 = rw*t2*t598;
        double t616 = rw*t3*t7*t598;
        double t618 = t190*t243*t549;
        double t623 = t2*t3*t321*t549;
        double t624 = t3*t7*t321*t549;
        double t627 = t190*t272*t549;
        double t630 = -t620;
        double t633 = t190*t243*(t276-t491)*(-1.0/2.0);

        double t636 = t3*t7*t321*(t276-t491)*(-1.0/2.0);
        double t637 = t361*(t276-t491)*(-1.0/2.0);
        double t643 = (t2*t3*t321*(t276-t491))/2.0;

        double t649 = t190*t272*(t276-t491)*(-1.0/2.0);
        double t650 = t8*t631;

        double t658 = t236*t598;
        double t662 = t3*t5*t113*t611;
        double t666 = (t2*t126*t611)/2.0;
        double t667 = (t7*t126*t611)/2.0;
        double t668 = (t399*t534)/2.0;
        double t672 = Itzz*t3*t74*t190*t598;
        double t678 = t399*(t276-t491)*(-1.0/4.0);
        double t679 = (t399*t548)/2.0;
        double t681 = (t236*t611)/2.0;
        double t682 = (t238*t611)/2.0;
        double t685 = (t3*t10*t190*t611)/2.0;
        double t689 = -rw*t2*(t446+(t2*t126*(t275-t415))/2.0);

        double t693 = t379+t599;
        double t698 = (t190*t243*t611)/2.0;
        double t708 = (t190*t272*t611)/2.0;

        double t716 = rw*t2*t700;
        double t719 = (t399*t598)/2.0;
        double t728 = ((t275-t415)*(t276-t491))/4.0;

        double t731 = (t399*t611)/4.0;
        double t733 = (t404*t611)/2.0;
        double t735 = -Iwyy*(t382-t402+(t5*(t246-t453))/2.0);

        double t747 = t526*(t276-t491)*(-1.0/2.0);

        double t756 = (t400*t680)/4.0;
        double t760 = (t598*(t275-t415))/2.0;
        double t767 = t323+t431+t604;
        double t813 = (t549*t611)/2.0;

        double t815 = rw*t2*t8*t802;
        double t816 = rw*t3*t7*t802;
        double t817 = rw*t7*t8*t802;
        double t820 = (t611*(t276-t491))/4.0;
        double t830 = (t526*t680)/4.0;
        double t834 = t117+t133+t206+t325+t343+t440+t558;

        double t843 = ((t441+t527)*(t276-t491))/2.0;
        double t849 = t238*t824;
        double t850 = (t549*t680)/4.0;
        double t852 = ((t276-t491)*(t446+(t2*t126*(t275-t415))/2.0))/2.0;
        double t864 = t539+t703;
        double t870 = t531+t714;
        double t896 = (t611*t663)/2.0;
        double t905 = Itzz*t3*t5*(t696-t10*t113*(t276-t491));
        double t916 = rw*t2*t8*t898;
        double t917 = rw*t3*t7*t898;

        double t929 = t680*(t446+(t2*t126*(t275-t415))/2.0)*(-1.0/4.0);
        double t935 = rw*t2*t8*t930;
        double t936 = rw*t3*t7*t930;
        double t971 = -Iwyy*t8*(-t575+t642+t190*t243*(t330-t416));

        double t976 = -Itzz*t3*t5*(-t575+t642+t190*t243*(t330-t416));

        double t979 = -rw*t3*t7*(-t575+t642+t190*t243*(t330-t416));
        double t980 = rw*t2*(-t575+t642+t190*t243*(t330-t416));
        double t982 = t655+t824;

        double t989 = t652+t828;
        double t990 = t671+t822;

        double t1086 = -rw*t2*(t832+t238*(t414-t744));
        double t1090 = Iwyy*t8*(t832+t238*(t414-t744));
        double t1092 = rw*t2*(t832+t238*(t414-t744));

        double t544 = t8*t536;
        double t579 = t52*t561;
        double t582 = t574*dom;

        double t588 = t283+t541;

        double t594 = Itzz*t10*t399*t476*(-1.0/2.0);
        double t595 = t191+t565;
        double t596 = t256+t550;
        double t597 = Iwyy*t590;
        double t609 = t8*t603;
        double t610 = t361+t536;

        double t613 = t2*t126*t590;
        double t614 = t7*t126*t590;
        double t615 = t8*t608;
        double t628 = -t618;
        double t629 = t336+t556;

        double t639 = rw*t2*t3*t607;
        double t640 = rw*t7*t8*t607;
        double t648 = rw*t7*(t331-t545);

        double t665 = t302*t590;
        double t673 = -t667;

        double t701 = Itzz*t10*t693;

        double t710 = -Iwyy*(t308-t497+rw*t3*t7*(t198-t363));
        double t722 = t8*t716;

        double t736 = -rw*t7*(t475-t567);
        double t742 = rw*t2*(t566-t238*(t198-t363));
        double t743 = rw*t7*t8*t87*(t475-t567);
        double t748 = t385+t682;
        double t751 = t473+t623;
        double t754 = t473+t633;
        double t757 = t444+t666;
        double t758 = t3*t74*t87*t190*(t475-t567);
        double t761 = t495+t608;

        double t765 = Itzz*t3*t74*t190*(t566-t238*(t198-t363));

        double t774 = t3*t767;
        double t775 = t489+t649;
        double t778 = t190*t243*(t566-t238*(t198-t363));

        double t797 = t236*(t489-t624);
        double t799 = t3*t10*t190*(t489-t624);
        double t800 = t559+t589;

        double t835 = t253+t333+t418+t428+t542+Tneta;
        double t841 = rw*t7*(t616-(t7*t126*(t276-t491))/2.0);
        double t846 = -t843;
        double t847 = (Iwyy*t834)/4.0;
        double t853 = (Itzz*t10*t834)/4.0;
        double t854 = ((t275-t415)*(t475-t567))/2.0;

        double t859 = (t2*t126*t834)/4.0;
        double t860 = (t7*t126*t834)/4.0;

        double t865 = ((t566-t238*(t198-t363))*(t275-t415))/2.0;
        double t873 = t8*t321*(t627-(t399*(t198-t363))/2.0);
        double t875 = rw*t7*t864;
        double t876 = t584+t678;
        double t879 = (t236*t834)/4.0;
        double t880 = (t238*t834)/4.0;
        double t887 = t432+t560+t568;
        double t888 = rw*t2*t8*t870;
        double t889 = rw*t3*t7*t870;

        double t903 = t238*(t685-(t8*t321*(t275-t415))/2.0);
        double t904 = t573+t708;
        double t906 = ((t275-t415)*(t489-t624))/2.0;
        double t907 = (t190*t243*t834)/4.0;
        double t910 = (t2*t3*t321*t834)/4.0;
        double t911 = (t3*t7*t321*t834)/4.0;
        double t913 = (t190*t272*t834)/4.0;

        double t942 = t658+t747;
        double t947 = t643+t763;
        double t952 = t555+t830;
        double t959 = -t7*t126*(t636+(t680*(t198-t363))/4.0);
        double t985 = t593+t850;
        double t993 = t681+t810;

        double t995 = rw*t7*t8*t982;
        double t996 = Iwyy*t990;
        double t997 = t97+t289+t312+t339+t445+t488+t513+t519;
        double t1000 = t3*t5*t87*t982;
        double t1004 = Itzz*t3*t5*t989;
        double t1011 = t399*(t636+(t680*(t198-t363))/4.0)*(-1.0/2.0);
        double t1015 = (t680*t834)/1.6e+1;
        double t1027 = -Itzz*t10*(-t654+t698+(t2*t3*t321*(t275-t415))/2.0);
        double t1030 = -rw*t2*t3*(-t654+t698+(t2*t3*t321*(t275-t415))/2.0);
        double t1031 = t543+t676+t731;
        double t1033 = ((t636+(t680*(t198-t363))/4.0)*(t275-t415))/2.0;
        double t1073 = t840+t849;
        double t1074 = t631+t980;
        double t1091 = t8*t1086;
        double t1104 = t263+t359+t365+t482+t486+t494+t502+t535+t600;
        double t587 = -t582;
        double t602 = t8*t597;

        double t625 = t2*t126*t595;
        double t638 = Iwyy*t629;
        double t644 = rw*t2*t3*t610;
        double t646 = rw*t7*t8*t610;

        double t661 = t10*t87*t629;
        double t664 = t7*t126*t629;
        double t669 = t2*t126*t629;
        double t697 = t190*t243*t629;
        double t706 = t190*t272*t629;
        double t711 = (t399*t588)/2.0;
        double t718 = (t399*t596)/2.0;
        double t745 = Itzz*t8*t742;
        double t766 = rw*t2*t757;
        double t768 = t2*t126*t751;

        double t771 = rw*t3*t7*t757;
        double t772 = rw*t7*t8*t761;
        double t776 = rw*t2*t3*t761;
        double t777 = Iwyy*t8*(t443+t673);
        double t781 = t236*t751;

        double t783 = t3*t10*t190*t751;
        double t784 = t238*t754;
        double t785 = rw*t2*t8*t775;
        double t788 = t236*t757;
        double t793 = t302*t751;
        double t798 = t190*t243*t748;

        double t804 = t190*t272*t748;
        double t806 = t238*t775;
        double t807 = t546+t628;
        double t809 = t190*t272*t751;
        double t812 = t190*t243*t757;
        double t818 = Itzz*t3*t5*t800;

        double t829 = t441+t527+t544;
        double t845 = (t399*t757)/2.0;

        double t866 = -t860;
        double t871 = t614+t630;
        double t881 = -t880;
        double t884 = t613+t657;

        double t899 = t2*t126*t876;
        double t900 = t7*t126*t876;
        double t902 = t578+t701;
        double t908 = -t903;

        double t912 = t3*t5*t10*t30*t887;
        double t914 = -t910;
        double t920 = rw*t3*t7*t904;
        double t922 = t302*t876;
        double t923 = t443+t538+t673;
        double t938 = t190*t243*(t347-t847);
        double t945 = rw*t2*t942;

        double t949 = Iwyy*t947;
        double t950 = rw*t2*(t859-t302*(t198-t363));
        double t953 = t2*t126*t947;
        double t955 = t328+t351+t514+t774;
        double t958 = t236*(t859-t302*(t198-t363));

        double t961 = rw*t2*t8*t952;
        double t962 = rw*t7*t8*t952;
        double t963 = t238*t947;
        double t965 = t3*t5*t87*t952;
        double t970 = t302*t947;
        double t975 = t190*t243*(t859-t302*(t198-t363));
        double t984 = t190*t272*t947;
        double t987 = Itzz*t985;
        double t988 = Iwyy*t985;
        double t1001 = t8*t996;
        double t1002 = t2*t126*t985;
        double t1003 = t7*t126*t985;
        double t1007 = (t399*(t859-t302*(t198-t363)))/2.0;
        double t1008 = (t399*t947)/2.0;
        double t1013 = t662+t853;
        double t1014 = t238*t993;

        double t1018 = t190*t243*t985;
        double t1019 = t190*t272*t985;
        double t1023 = t656+t904;
        double t1024 = t622+t913;
        double t1036 = Iwyy*t1031;
        double t1042 = t728+t879;
        double t1043 = t985*(t275-t415)*(-1.0/2.0);
        double t1058 = -t7*t126*(t911+(t611*(t198-t363))/2.0);
        double t1064 = -t236*(t911+(t611*(t198-t363))/2.0);
        double t1066 = t799+t873;
        double t1071 = t813+t879;
        double t1076 = t663+t736+t742;
        double t1078 = Iwyy*t8*t1073;
        double t1083 = rw*t2*t3*t1074;
        double t1084 = rw*t7*t8*t1074;
        double t1095 = t399*(t911+(t611*(t198-t363))/2.0)*(-1.0/2.0);
        double t1101 = -rw*t7*(-t722+t907+(t404*(t275-t415))/2.0);
        double t1106 = -rw*t2*(-t719+t778+t190*t272*(t475-t567));
        double t1120 = t820+t1015;
        double t1170 = t1000+t1027;
        double t1173 = t639+t713+t827+t875;
        double t1182 = t979+t1073;

        double t1204 = t10*t87*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567)));
        double t1206 = -rw*t7*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567)));
        double t659 = -t644;
        double t674 = -t664;
        double t675 = t37+t300+t587;
        double t702 = t345+t638;
        double t707 = -t697;
        double t838 = rw*t7*t829;
        double t848 = -t845;
        double t862 = t8*t321*t807;
        double t863 = t583+t669;
        double t882 = Iwyy*t8*t871;
        double t883 = rw*t7*t871;
        double t892 = Itzz*t3*t5*t871;
        double t894 = Iwyy*t8*t884;
        double t895 = rw*t2*t884;
        double t901 = Itzz*t3*t5*t884;

        double t924 = t238*t902;
        double t925 = rw*t7*t923;
        double t933 = t679+t706;
        double t943 = t487+t866;
        double t948 = t8*t945;
        double t957 = t955*dth;
        double t991 = t672+t818;

        double t1006 = -t1003;
        double t1009 = -t1008;
        double t1010 = t601+t881;
        double t1016 = Itzz*t10*t1013;
        double t1020 = t756+t771;
        double t1028 = rw*t2*t1024;
        double t1029 = Itzz*t10*t1023;
        double t1034 = rw*t3*t7*t1024;
        double t1035 = t665+t908;
        double t1040 = t8*t1036;
        double t1044 = t238*t1042;
        double t1046 = t733+t914;
        double t1067 = Iwyy*t1066;
        double t1068 = rw*t2*t1066;

        double t1080 = t2*t126*t1071;
        double t1081 = t7*t126*t1071;
        double t1082 = -rw*t2*(t806-t899);
        double t1093 = t236*t1071;
        double t1094 = t661+t987;
        double t1096 = t190*t243*t1071;
        double t1097 = t190*t272*t1071;
        double t1107 = t8*t1106;

        double t1110 = t852+t958;
        double t1111 = -rw*t2*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624));
        double t1113 = -rw*t3*t7*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624));
        double t1117 = -t3*t10*t190*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624));
        double t1125 = t87*t1120;
        double t1126 = -Itzz*t10*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673));
        double t1127 = -Iwyy*t8*(t798+(t399*(t443+t673))/2.0+(t498*(t275-t415))/2.0);
        double t1130 = -rw*t2*t3*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673));
        double t1134 = t2*t126*t1120;
        double t1135 = t7*t126*t1120;
        double t1141 = t190*t243*t1120;
        double t1142 = t190*t272*t1120;
        double t1146 = ((t275-t415)*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624)))/2.0;
        double t1148 = (t399*t1120)/2.0;
        double t1152 = -rw*t3*t7*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624));
        double t1156 = t953+t959;
        double t1159 = ((t276-t491)*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673)))/2.0;

        double t1163 = Itzz*(t1002-t238*(t636+(t680*(t198-t363))/4.0));
        double t1165 = rw*t2*(t1002-t238*(t636+(t680*(t198-t363))/4.0));
        double t1166 = t625+t710+t776+t841;
        double t1168 = t533+t637+t816+t945;
        double t1180 = rw*t7*t8*(t650+rw*t2*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)));
        double t1183 = rw*t7*t1182;
        double t1196 = -rw*t7*t8*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0));
        double t1211 = -rw*t7*(t798+(t399*(t443+t673))/2.0+t8*t980+(t498*(t275-t415))/2.0);
        double t1213 = -rw*t2*t3*(t798+(t399*(t443+t673))/2.0+t8*t980+(t498*(t275-t415))/2.0);
        double t1219 = t797+t1011+t1019;
        double t1233 = -rw*t2*t3*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)));
        double t1244 = -rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)));
        double t1291 = -Iwyy*(t1120+rw*t2*t8*(t636+(t680*(t198-t363))/4.0)+rw*t3*t7*(t911+(t611*(t198-t363))/2.0));
        double t1300 = -rw*t7*(t961+t236*(t443+t673)+(t680*(t441+t527))/4.0+rw*t3*t7*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673)));
        double t1323 = t842+t920+t993+t995+t1030;

        double t686 = t675*dph*2.0;
        double t737 = t3*t10*t190*t702;
        double t750 = t190*t243*t702;
        double t752 = t190*t272*t702;
        double t844 = -t838;
        double t855 = t571+t674;
        double t872 = rw*t2*t863;
        double t886 = t10*t87*t863;
        double t928 = -t925;
        double t931 = t668+t707;
        double t940 = rw*t2*t933;
        double t954 = t236*t943;
        double t973 = t190*t272*t943;
        double t981 = t358+t469+t648+t659;
        double t999 = (t399*t943)/2.0;
        double t1012 = Itzz*t1010;
        double t1022 = Iwyy*t8*t1020;
        double t1025 = t594+t924;
        double t1032 = t190*t243*t1010;
        double t1038 = t190*t272*t1010;
        double t1045 = t615+t943;
        double t1047 = -t1044;
        double t1048 = Iwyy*t1046;
        double t1052 = t758+t892;
        double t1053 = t2*t126*t1046;

        double t1057 = t765+t901;
        double t1060 = t236*t1046;
        double t1061 = t238*t1046;
        double t1069 = t190*t272*t1046;
        double t1070 = -t1068;
        double t1087 = -t1081;
        double t1088 = (t399*t1046)/2.0;
        double t1089 = t8*t1082;
        double t1099 = t302*t1094;
        double t1112 = t8*t1111;
        double t1114 = rw*t2*t1110;
        double t1116 = Itzz*t1113;
        double t1123 = t694+t804+t848;
        double t1138 = -t1135;
        double t1139 = t597+t883+t895;
        double t1149 = -t1148;

        double t1157 = rw*t2*t1156;
        double t1160 = t963+t1006;
        double t1167 = Itzz*t3*t5*t1166;
        double t1171 = rw*t7*t8*t1168;
        double t1172 = t1004+t1029;
        double t1184 = t533+t720+t1111;
        double t1185 = t8*t1183;
        double t1186 = -t1183;
        double t1208 = t965+t1126;
        double t1210 = t781+t1009+t1018;
        double t1216 = rw*t2*(t1134-t302*(t636+(t680*(t198-t363))/4.0));
        double t1221 = Iwyy*t1219;
        double t1222 = rw*t2*t1219;
        double t1225 = t7*t126*t1219;
        double t1227 = t302*t1219;
        double t1240 = t785+t1034+t1042;
        double t1246 = t8*t1244;
        double t1269 = t190*t243*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624));

        double t1271 = -rw*t7*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567)));
        double t1272 = t905+t1016+t1125;
        double t1275 = t906+t1095+t1097;
        double t1294 = -rw*t3*t7*(Itzz*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624))+t3*t5*t87*t863);
        double t1297 = t788+t929+t962+t1130;
        double t1303 = t1033+t1064+t1142;
        double t1315 = -t7*t126*(-t1093+t1148+(t985*(t275-t415))/2.0);
        double t1316 = t7*t126*(-t1093+t1148+(t985*(t275-t415))/2.0);
        double t1324 = Iwyy*t1323;

        double t868 = Itzz*t10*t855;
        double t869 = rw*t7*t855;
        double t877 = Itzz*t3*t5*t855;
        double t932 = rw*t7*t931;
        double t941 = -t940;
        double t998 = t3*t5*t87*t981;
        double t1039 = -t1032;
        double t1041 = t3*t5*t87*t1025;
        double t1049 = rw*t7*t1045;
        double t1055 = rw*t7*t1052;
        double t1063 = rw*t2*t1057;
        double t1115 = Itzz*t1112;
        double t1128 = rw*t2*t1123;
        double t1129 = Iwyy*t8*t1123;
        double t1137 = rw*t3*t7*t1123;
        double t1140 = t603+t646+t689+t844;
        double t1158 = t8*t1157;

        double t1164 = t87*t1160;
        double t1169 = t350+t743+t745+t1012;
        double t1175 = t640+t735+t766+t928;
        double t1188 = rw*t2*t3*t1184;
        double t1189 = rw*t7*t8*t1184;
        double t1194 = t760+t973+t975;
        double t1212 = Iwyy*t1210;
        double t1215 = t886+t1163;
        double t1217 = t2*t126*t1210;
        double t1218 = -t1216;
        double t1223 = t8*t1222;
        double t1224 = t302*t1210;
        double t1237 = t865+t1007+t1038;
        double t1241 = Iwyy*t1240;
        double t1247 = t749+t1053+t1058;
        double t1253 = t793+t1061+t1087;
        double t1268 = rw*t2*t3*(t836-t1069+t190*t243*(t911+(t611*(t198-t363))/2.0));

        double t1276 = Iwyy*t1275;
        double t1277 = t238*t1272;
        double t1278 = rw*t2*t1275;
        double t1279 = rw*t3*t7*t1275;
        double t1282 = t718+t752+t1206;
        double t1285 = -rw*t7*(t1160+rw*t3*t7*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624)));
        double t1296 = rw*t7*t8*(t712+t889-t949-t1157);
        double t1299 = rw*t2*t1297;
        double t1304 = -t238*(-t1060+t1141+(t947*(t275-t415))/2.0);
        double t1305 = t238*t1303;
        double t1311 = t1043+t1093+t1149;
        double t1321 = t1152+t1210;
        double t1336 = -rw*t7*(-t1088+t1096+(t751*(t275-t415))/2.0+rw*t2*t8*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)));
        double t1337 = t894+t1067+t1271;
        double t1344 = t996+t1083+t1092+t1186;
        double t874 = -t869;
        double t1133 = -t1128;
        double t1144 = Itzz*t10*t1140;
        double t1154 = t8*t321*t1140;
        double t1155 = t450+t932+t941;
        double t1181 = t3*t10*t190*t1175;

        double t1197 = rw*t2*t1194;
        double t1198 = rw*t7*t1194;
        double t1220 = rw*t2*t8*t1215;
        double t1236 = t854+t999+t1039;
        double t1238 = rw*t2*t1237;
        double t1239 = rw*t3*t7*t1237;
        double t1242 = (t680*t1194)/4.0;
        double t1248 = rw*t2*t1247;
        double t1249 = rw*t3*t7*t1247;
        double t1250 = t236*t1247;
        double t1252 = (t399*t1247)/2.0;
        double t1254 = Itzz*t1253;
        double t1258 = t522+t737+t1055+t1063;
        double t1262 = t190*t272*t1253;
        double t1280 = -t1278;
        double t1287 = -Itzz*(t347+t529-t772-t847-t950+t1049);
        double t1306 = t868+t1116+t1164;
        double t1322 = rw*t7*t1321;
        double t1325 = t867+t1014+t1091+t1137;
        double t1327 = t1112+t1253;
        double t1345 = Itzz*t3*t5*t1344;
        double t1352 = -Itzz*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)));
        double t1356 = rw*t7*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)));
        double t1362 = t1196+t1268+t1303;
        double t1364 = t3*t10*t190*(-t725+t917+t988+t1165+t1188+t1285);
        double t1368 = t1223+t1279+t1311;
        double t1382 = t1299+t1300+t1324;
        double t1153 = t702+t872+t874;
        double t1199 = t3*t1198;
        double t1200 = -t1197;

        double t1230 = t998+t1144;
        double t1251 = -t1250;
        double t1260 = t3*t5*t87*t1258;
        double t1263 = -t1262;
        double t1288 = t10*t1287;
        double t1301 = t1154+t1181;
        double t1307 = rw*t7*t8*t1306;
        double t1317 = t1107+t1236;
        double t1326 = Iwyy*t8*t1325;
        double t1329 = rw*t7*t1327;
        double t1338 = t877+t1115+t1254;
        double t1339 = t777+t888+t1048+t1248;
        double t1346 = t1036+t1084+t1133+t1211;
        double t1348 = t970+t1138+t1158+t1249;
        double t1360 = t1185+t1213+t1325;
        double t1367 = t1204+t1352;
        double t1369 = Iwyy*t1368;
        double t1381 = t1040+t1180+t1280+t1336;
        double t1394 = t1090+t1221+t1233+t1356;
        double t1232 = Itzz*t3*t5*t1230;
        double t1313 = rw*t2*t3*(t815+t938+t1200-(t495*(t275-t415))/2.0);
        double t1318 = rw*t7*t1317;
        double t1319 = rw*t2*t3*t1317;

        double t1331 = t846+t948+t954+t1199;
        double t1340 = rw*t2*t3*t1338;
        double t1341 = rw*t2*t3*t1339;
        double t1342 = t1167+t1288;
        double t1347 = -t10*t87*t1346;
        double t1349 = rw*t7*t1348;
        double t1350 = t1159+t1242+t1251;
        double t1372 = t1146+t1252+t1263+t1269;
        double t1384 = -t3*t10*t190*(-t896+t916-t1189+t1329+t236*(t347-t847)+rw*t2*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)));

        double t1332 = rw*t7*t1331;
        double t1374 = rw*t2*t1372;
        double t1375 = rw*t7*t1372;
        double t1386 = t1345+t1347;
        double t1389 = t922+t1047+t1089+t1239+t1246+t1319;
        double t1405 = -Itzz*(t1384+t8*t321*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))));
        double t1406 = t1022+t1218+t1291+t1296+t1341+t1349;
        double t1333 = -t1332;
        double t1376 = t3*t1375;
        double t1390 = Itzz*t1389;
        double t1408 = t1260+t1405;
        double t1396 = Itzz*(t1114+t1171+t1241+t1313+t1333-Iwyy*t8*(t467+rw*t3*t7*(t446+(t2*t126*(t275-t415))/2.0)));
        double t1409 = rw*t8*t1408;
        double t1414 = -rw*t7*(-t1224+t1315+t1376+t238*(-t1060+t1141+(t947*(t275-t415))/2.0)+rw*t2*t8*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0))));
        double t1410 = t1409*4.0;

        double et1 = t1409;
        double et2 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et3 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double t1421 = -1.0/(et1+et2+et3);
        double et4 = t1410;
        double et5 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et6 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double t1422 = -1.0/(et4+et5+et6);

        double et7 = t1410;
        double et8 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et9 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et10 = t1409;
        double et11 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et12 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et13 = t1409;
        double et14 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et15 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et16 = t1409;
        double et17 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et18 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et19 = t1422*(Itzz*t1344+t30*t74*t887+Itzz*rw*t8*t1139)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352)))+t835*t1421*(t1345-Itzz*t10*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+rw*t8*t1258);
        double et20 = t1104*t1422*(t912+Itzz*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))-Itzz*rw*t3*t8*t10*t190*t1076)+(t997*(Itzz*(-t725+t917+t988+t1165+t1188+t1285)+t10*t87*t1153+Itzz*rw*t72*t321*t1076))/(et7+et8+et9);
        double et21 = (t626*(Itzz*(t1001-t1222-t1322+rw*t2*t3*(t650+rw*t2*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))-Itzz*t10*t1155+rw*t8*t87*(-t602+t1068+rw*t7*(t783-t862))))/(et10+et11+et12);
        double et22 = (rw*t368*(-Itzz*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+t10*t87*(-t711+t750+rw*t2*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567))))+Itzz*rw*t8*(t882+Iwyy*(t783-t862)+rw*t2*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))))))/(et13+et14+et15)+rw*t376*t1421*(Itzz*t1394+t10*t87*t1282+Itzz*rw*t8*t1337);
        double et23 = (rw*t3*t52*(Itzz*(t1364+t8*t321*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t10*t1258))/(et16+et17+et18);
        double et24 = t1410;
        double et25 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et26 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et27 = t1409;
        double et28 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et29 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et30 = t1409;
        double et31 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et32 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et33 = t1409;
        double et34 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et35 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et36 = t1409;
        double et37 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et38 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et39 = ((t912-Itzz*t1346)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et24+et25+et26)+(t835*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+t3*t5*t87*t1346))/(et27+et28+et29)+t1104*t1422*(Itzz*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))-t30*t59*t61*t887)+(t626*(Itzz*t1381+Itzz*t3*t5*t1155))/(et30+et31+et32);
        double et40 = t997*t1422*(Itzz*(-t896+t916-t1189+t1329+t236*(t347-t847)+rw*t2*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)))+t3*t5*t87*t1153)+(rw*t368*(Itzz*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))-t3*t5*t87*(-t711+t750+rw*t2*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567))))))/(et33+et34+et35);
        double et41 = (rw*t376*(Itzz*(t1129-t1276-t1375+rw*t7*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+t3*t5*t87*t1282))/(et36+et37+et38)+rw*t3*t52*t1408*t1421;
        double et42 = t1410;
        double et43 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et44 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et45 = t1409;
        double et46 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et47 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et48 = t1409;
        double et49 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et50 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et51 = t87*(-t1093+t1148+(t985*(t275-t415))/2.0)+Itzz*t10*((t399*t1013)/2.0+(t629*(t275-t415))/2.0)-rw*t8*(Itzz*((t590*t834)/4.0-t549*(t685-(t8*t321*(t275-t415))/2.0))+Itzz*t8*t1070-Itzz*t3*t5*(t3*t10*t190*t629+t3*t5*t113*t321*t415)+rw*t7*t8*t87*(t783-t862))-rw*t2*t3*(Itzz*(-t1088+t1096+(t751*(t275-t415))/2.0)-Itzz*t3*t5*t931+Itzz*rw*t2*t8*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))+rw*t3*t7*(Itzz*t1275-t3*t5*t87*t933);
        double et52 = Itzz*t3*t5*(t236*t629+(t399*(t696-t10*t113*(t276-t491)))/2.0)+rw*t7*t8*(Itzz*t1152+Itzz*t1210+Itzz*t10*t931)+rw*t2*t8*(Itzz*t1219+t10*t87*t933);
        double et53 = t997*t1422*(t1099+t1220+t1277+t1294+t1307+t1340+rw*t72*t321*t1169)+((t87*t1360+Itzz*t10*t1025+rw*t8*(Itzz*t1035+Itzz*t8*t883+Itzz*t8*t895-rw*t30*t59*t61*t74*t204))*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et42+et43+et44)+t1104*t1422*(t1041+t1390+rw*t3*t8*t10*t190*t1169)+(t835*(t10*t1390+rw*t8*(t8*t1055+t8*t1063+Itzz*t3*t5*t1035+t3*t74*t190*t1012)+t3*t5*t87*t1360))/(et45+et46+et47)+(t626*(et51+et52))/(et48+et49+et50);
        double et54 = rw*t376*t1421*(t87*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+Itzz*t10*(Itzz*t10*t1237+t3*t5*t87*t1123)+rw*t8*(t87*(t8*t321*t1237+t3*t10*t190*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)))+Itzz*t3*t5*t1057+Itzz*rw*t7*t8*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))))+rw*t7*t8*t1367+rw*t2*t3*(Itzz*t1372-Itzz*t3*t5*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567))))-Itzz*t3*t5*(Itzz*t10*(t806-t899)+t3*t5*t87*(t832+t238*(t414-t744))));
        double et55 = rw*t368*t1421*(Itzz*(t1224+t1304+t1316)-t10*t87*(Itzz*t10*t1236+Itzz*t3*t5*(t798+(t399*(t443+t673))/2.0+(t498*(t275-t415))/2.0))+rw*t8*(-Itzz*(t8*t321*t1236-t3*t10*t190*t1253)+t3*t5*t87*t1052+Itzz*rw*t2*t8*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))))-t3*t5*t87*(Itzz*t10*(t784-t900)-t3*t5*t87*t1073)+rw*t2*t8*t1367-rw*t3*t7*(Itzz*t1372-Itzz*t3*t5*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567)))))+rw*t3*t122*t1421*(t8*t321*(t1041+t1390)-t3*t10*t190*(t1099+t1220+t1277+t1294+t1307+t1340));
        double et56 = t1410;
        double et57 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et58 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et59 = -Itzz*(-t8*t1324+rw*t2*t1362+rw*t7*(-t1060+t1141+(t947*(t275-t415))/2.0+rw*t2*t8*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0))+rw*t3*t7*(t836-t1069+t190*t243*(t911+(t611*(t198-t363))/2.0))));
        double et60 = rw*t8*(Itzz*(t8*t321*(t609-t1028+t1101+rw*t7*t8*(t369-t716))-t3*t10*t190*(Iwyy*t8*(t382-t402+(t5*(t246-t453))/2.0)-rw*t2*(t911+(t611*(t198-t363))/2.0)+rw*t7*(-t733+t910+rw*t2*t8*(t499+t2*t3*t321*(t198-t363)))-rw*t7*t8*(t8*t362+rw*t2*(t499+t2*t3*t321*(t198-t363)))))-t3*t5*t87*(rw*t7*(t3*t10*t190*t534+Itzz*t3*t8*t228*t243*t321)-rw*t2*(t3*t10*t190*t548+t3*t8*t87*t228*t272*t321)+Itzz*Iwyy*t3*t8*t74*t190));
        double et61 = t10*t87*(Itzz*t10*(t609-t1028+t1101+rw*t7*t8*(t369-t716))-Itzz*t3*t5*(t8*t469+rw*t7*(t754+rw*t3*t7*t700)+rw*t2*t775-rw*t2*t3*(t369-t716)))+Itzz*t3*t5*(rw*t2*t1172-rw*t7*t1170);
        double et62 = t1409;
        double et63 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et64 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et65 = t1410;
        double et66 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et67 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et68 = Itzz*(Iwyy*(-t1060+t1141+(t947*(t275-t415))/2.0+rw*t2*t8*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0))+rw*t3*t7*(t836-t1069+t190*t243*(t911+(t611*(t198-t363))/2.0)))+rw*t2*t1350-Iwyy*t8*(t961+t236*(t443+t673)+(t680*(t441+t527))/4.0+rw*t3*t7*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673))))-rw*t8*(Itzz*(t8*t321*(t815+t938+t1200-(t495*(t275-t415))/2.0)+t3*t10*t190*t1339)+Itzz*t3*t5*(rw*t2*t991+t3*t10*t190*t588+Itzz*Iwyy*t3*t8*t228*t243*t321));
        double et69 = Itzz*t10*(Itzz*t10*(t815+t938+t1200-(t495*(t275-t415))/2.0)+t3*t5*t87*t1168)+Itzz*t3*t5*(Iwyy*t1170+rw*t2*t1208);
        double et70 = t1409;
        double et71 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et72 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et73 = rw*t376*(Itzz*(Iwyy*t1362+Iwyy*t8*t1297-rw*t7*t1350)-t10*t87*(Itzz*t10*(-t817+t1198+(t516*(t275-t415))/2.0+t190*t272*(t347-t847))-t3*t5*t87*(t551+rw*t7*t942-rw*t2*t3*t802+(Iwyy*t190*t272*(t276-t491))/2.0))-rw*t8*(Itzz*(t8*t321*(-t817+t1198+(t516*(t275-t415))/2.0+t190*t272*(t347-t847))-t3*t10*t190*(-Iwyy*(t911+(t611*(t198-t363))/2.0)+Iwyy*t8*t757+rw*t7*t1247+rw*t7*t8*t870))-Itzz*t3*t5*(rw*t7*t991+t3*t10*t190*t596+Iwyy*t3*t8*t87*t228*t272*t321))+t3*t5*t87*(Iwyy*t1172+rw*t7*t1208));
        double et74 = 1.0/(et70+et71+et72);
        double et75 = (t997*(Itzz*t1406-t10*t87*t1342+rw*t72*t321*(t426+Itzz*(t347+t529-t772-t847-t950+t1049))+Itzz*t3*t5*(Itzz*t10*t1175-t3*t5*t87*t1173)))/(et56+et57+et58)+t626*t1421*(et59+et60+et61)+t1104*t1422*(t1232+t1396-rw*t3*t8*t10*t190*(t426+Itzz*(t347+t529-t772-t847-t950+t1049)))+(t835*(t10*t1396+rw*t8*(Itzz*t3*t5*t1301+t3*t74*t190*t1287)+Itzz*t3*t5*t1382))/(et62+et63+et64)+((Itzz*t1382+rw*t8*(Itzz*t1301+t30*t59*t74*t228*t370)+t10*t87*t1230)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et65+et66+et67)+rw*t368*t1421*(et68+et69);
        double et76 = et73*et74+rw*t3*t122*t1421*(t8*t321*(t1232+t1396)+t3*t10*t190*(Itzz*t1406-t10*t87*t1342+Itzz*t3*t5*(Itzz*t10*t1175-t3*t5*t87*t1173)));
        double et77 = t1410;
        double et78 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et79 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);

        double et80 = t1410;
        double et81 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et82 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et83 = t1409;
        double et84 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et85 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et86 = t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))-rw*t8*(t1384+t8*t321*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))));
        double et87 = rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))));
        double et88 = t1409;
        double et89 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et90 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et91 = t10*t87*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+Itzz*t3*t5*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))));
        double et92 = -Itzz*rw*t3*t5*t8*(t882+Iwyy*(t783-t862)+rw*t2*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))));
        double et93 = t1409;
        double et94 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et95 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et96 = t1409;
        double et97 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et98 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et99 = t1409;
        double et100 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et101 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et102 = ((t1386+Itzz*rw*t3*t5*t8*t1139)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et77+et78+et79);
        double et103 = (t1104*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))-Itzz*t8*t10*t59*t238*t1076))/(et80+et81+et82)+(t835*(et86+et87))/(et83+et84+et85);
        double et104 = (t626*(-Itzz*t10*t1381+t3*t5*t87*(t1001-t1222-t1322+rw*t2*t3*(t650+rw*t2*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+Itzz*rw*t3*t5*t8*(-t602+t1068+rw*t7*(t783-t862))))/(et88+et89+et90)+t997*t1422*(t10*t87*(-t896+t916-t1189+t1329+t236*(t347-t847)+rw*t2*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)))+Itzz*t3*t5*(-t725+t917+t988+t1165+t1188+t1285)+Itzz*rw*t3*t5*t72*t321*t1076)+(rw*t368*(et91+et92))/(et93+et94+et95);
        double et105 = (rw*t376*(t10*t87*(t1129-t1276-t1375+rw*t7*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+Itzz*t3*t5*t1394+Itzz*rw*t3*t5*t8*t1337))/(et96+et97+et98);
        double et106 = (rw*t3*t122*(Itzz*t10*(t1384+t8*t321*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))+Itzz*t3*t5*(t1364+t8*t321*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))))/(et99+et100+et101);

        d_dPS[index] = et53+et54+et55;
        
        index = index + num_threads;
        // index1 = (grid_size[0] > 1) ? index : 0;
        index2 = (grid_size[1] > 1) ? index : 0;
        index3 = (grid_size[2] > 1) ? index : 0;
        index4 = (grid_size[3] > 1) ? index : 0;
        index5 = (grid_size[4] > 1) ? index : 0;
        index6 = (grid_size[5] > 1) ? index : 0;
        index7 = (grid_size[6] > 1) ? index : 0;
        index8 = (grid_size[7] > 1) ? index : 0;

        uindex1 = (active_actions[0] == 1) ? index : 0;
        uindex2 = (active_actions[1] == 1) ? index : 0;
    }
}

__global__ void dyn6_mex_continuous(double* const d_dOM, // outputs
    double* const PH, double* const TH, double* const OM, double* const dPH, double* const dTH,
    double* const dPS, double* const dOM, double* const dNE, // input states
    double const* const in1, double const* const in2, // input actions
    const double mw, const double Iwxx, const double Iwyy, const double Iwzz, 
    const double mf, const double Ifxx, const double Ifyy, const double Ifzz, 
    const double mt, const double Itxx, const double Ityy, const double Itzz, 
    const double nw, const double nt, const double rw, const double rf, const double rt, const double fcoeff, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3]
                        * grid_size[4] * grid_size[5] * grid_size[6] * grid_size[7];
    double th, ps, om, ne, dph, dth, dps, dom, dne, Tneta, Tomega;
    ps = 0; ne = 0;
    
    // int index1 = (grid_size[0] > 1) ? index : 0;
    int index2 = (grid_size[1] > 1) ? index : 0;
    int index3 = (grid_size[2] > 1) ? index : 0;
    int index4 = (grid_size[3] > 1) ? index : 0;
    int index5 = (grid_size[4] > 1) ? index : 0;
    int index6 = (grid_size[5] > 1) ? index : 0;
    int index7 = (grid_size[6] > 1) ? index : 0;
    int index8 = (grid_size[7] > 1) ? index : 0;

    int uindex1 = (active_actions[0] == 1) ? index : 0;
    int uindex2 = (active_actions[1] == 1) ? index : 0;

    while (index < num_elements)
    {
        th = TH[index2];
        om = OM[index3];
        dph = dPH[index4];
        dth = dTH[index5];
        dps = dPS[index6];
        dom = dOM[index7];
        dne = dNE[index8];

        Tomega = nw*in1[uindex1];
        Tneta = nt*in2[uindex2];

        double t2 = cos(ph);
        double t3 = cos(th);
        double t4 = cos(ps);
        double t5 = cos(om);
        double t6 = cos(ne);
        double t7 = sin(ph);
        double t8 = sin(th);
        double t9 = sin(ps);
        double t10 = sin(om);
        double t11 = sin(ne);
        double t12 = mf+mt;
        double t13 = mf*rf;
        double t14 = mt*rt;
        double t15 = rf+rt;
        double t16 = Ifxx*2.0;
        double t17 = Ifxx*4.0;
        double t18 = Ifyy*2.0;
        double t19 = Ifyy*4.0;
        double t20 = Ifzz*2.0;
        double t21 = Ifzz*4.0;
        double t22 = Itxx*2.0;
        double t23 = Itxx*3.0;
        double t24 = Itxx*4.0;
        double t25 = Ityy*2.0;
        double t26 = Ityy*3.0;
        double t27 = Ityy*4.0;
        double t28 = Itzz*2.0;
        double t29 = Itzz*4.0;
        double t30 = Itzz*Itzz;
        double t31 = Iwxx*2.0;
        double t32 = Iwxx*4.0;
        double t33 = Iwyy*2.0;
        double t34 = Iwyy*4.0;
        double t35 = Iwzz*2.0;
        double t36 = Iwzz*4.0;
        double t37 = fcoeff*2.0;

        double t40 = mt*2.0;
        double t42 = rf*rf;

        double t44 = rw*rw;
        double t45 = Tomega*4.0;
        double t46 = th*2.0;
        double t47 = ps*2.0;
        double t48 = om*2.0;
        double t49 = ne*2.0;
        double t50 = dom*2.0;
        double t51 = dph*dph;
        double t52 = dth*dth;
        double t53 = dom*dom;
        double t79 = Ifyy*8.0;
        double t80 = -Ifzz;
        double t83 = -Ityy;
        double t87 = -Itzz;
        double t92 = Iwyy*8.0;
        double t93 = -Iwzz;
        double t109 = Ifxx/2.0;
        double t110 = Ifzz/2.0;
        double t111 = Itxx/4.0;
        double t112 = Ityy/4.0;
        double t113 = Itzz/2.0;
        double t114 = Iwxx/2.0;
        double t115 = Iwzz/2.0;
        double t54 = cos(t46);
        double t55 = cos(t47);
        double t56 = cos(t48);
        double t57 = cos(t49);
        double t58 = t2*t2;
        double t59 = t3*t3;
        double t60 = t4*t4;
        double t61 = t5*t5;
        double t62 = t5*t5*t5;
        double t63 = t6*t6;
        double t64 = t13*2.0;
        double t65 = t13*4.0;
        double t66 = t14*2.0;
        double t67 = sin(t46);
        double t68 = sin(t47);
        double t69 = sin(t48);
        double t70 = sin(t49);
        double t71 = t7*t7;
        double t72 = t8*t8;
        double t73 = t9*t9;
        double t74 = t10*t10;
        double t75 = t11*t11;
        double t76 = mw+t12;
        double t77 = -t16;
        double t78 = -t17;
        double t81 = -t20;
        double t82 = -t21;
        double t84 = -t25;
        double t85 = -t26;
        double t86 = -t27;
        double t88 = -t28;
        double t89 = -t29;
        double t90 = -t31;
        double t91 = -t32;
        double t94 = -t35;
        double t95 = -t36;
        double t96 = t8*dph;
        double t97 = -t45;
        double t98 = rf*t12;
        double t99 = mt*t15;
        double t100 = t2*t5;
        double t101 = t2*t10;
        double t102 = t5*t7;
        double t103 = t6*t8;
        double t104 = t2*dph*dps;
        double t105 = t7*t10;
        double t106 = t8*t11;
        double t107 = rf*t13;
        double t108 = t15*t15;
        double t116 = rf*t14*4.0;
        double t117 = rf*t14*6.0;
        double t118 = Ifyy*Iwyy*t8;
        double t120 = t15*t40;
        double t122 = -t52;
        double t124 = Itxx+t83;
        double t125 = Iwxx+t93;
        double t133 = rt*t14*3.0;
        double t134 = rt*t14*4.0;
        double t135 = rf*t14*8.0;
        double t140 = t3*t6*t10;
        double t141 = t51+t52;
        double t144 = t3*t10*t11;
        double t148 = t20+t28;
        double t155 = t12*2.0;
        double t156 = t12*3.0;
        double t157 = t2*t3*dph*dth*2.0;
        double t158 = t7*t8*t52;
        double t159 = t3*t7*dph*dth*2.0;
        double t179 = Iwyy*mt*t8*t42;
        double t180 = Iwyy*rt*t8*t14;
        double t181 = (Itxx*Iwyy*t8)/2.0;
        double t182 = (Ityy*Iwyy*t8)/2.0;
        double t183 = rf*t8*t14*t33;
        double t191 = Iwyy*rw*t3*t7*t8;
        double t219 = t12*t42*4.0;
        double t119 = t98*2.0;
        double t121 = t99*4.0;
        double t123 = t50*t96;
        double t126 = rw*t76;
        double t127 = t96+dps;
        double t128 = t96+dom;
        double t129 = t54*3.0;
        double t130 = rf*t64;
        double t131 = rf*t65;
        double t132 = rt*t66;

        double t137 = Itxx*t63;
        double t138 = Ityy*t75;
        double t139 = t8*t100;
        double t142 = t8*t101;
        double t143 = t8*t102;
        double t145 = t8*t105;
        double t146 = t15*t99;
        double t149 = t22*t63;
        double t150 = t27*t63;
        double t151 = t35*t60;
        double t152 = t24*t75;
        double t153 = t25*t75;
        double t154 = t31*t73;
        double t160 = Iwyy*t13*t101;
        double t161 = Iwyy*mt*rf*t101;
        double t162 = Iwyy*t14*t101;
        double t163 = Iwyy*t13*t105;
        double t164 = Iwyy*mt*rf*t105;
        double t165 = Iwyy*t14*t105;
        double t166 = t40*t108;
        double t168 = t51*t54;
        double t169 = t51*t59;
        double t171 = t107/2.0;
        double t172 = t22+t84;
        double t173 = t23+t85;
        double t174 = t24+t86;
        double t175 = t124*t124;
        double t176 = t31+t94;
        double t177 = t32+t95;
        double t178 = Iwyy*t8*t107;
        double t184 = t50+t96;
        double t186 = -t140;
        double t187 = -t157;
        double t189 = t14+t98;
        double t190 = t13+t99;
        double t203 = t56*t59*2.0;
        double t205 = t42*t155;
        double t206 = t42*t156;

        double t208 = Ityy*Iwyy*t11*t140;
        double t210 = t64+t120;
        double t212 = t57*t124;
        double t213 = t55*t125;
        double t214 = t68*t125;
        double t215 = t44*t58*t76;

        double t217 = t57*t182;
        double t218 = t44*t71*t76;
        double t223 = t2*t8*t141;
        double t225 = t61*t148;
        double t226 = t10*t57*t67*4.0;
        double t232 = Itxx*Iwyy*t8*t57*(-1.0/2.0);
        double t236 = t5*t6*t11*t124;
        double t244 = t103+t144;
        double t246 = t8*t70*t124;
        double t253 = t6*t11*t53*t124;
        double t254 = t4*t9*t52*t125;

        double t275 = t3*t10*t70*t124;
        double t287 = t5*t67*t70*t124;
        double t288 = t3*t69*t70*t124;

        double t333 = t6*t11*t61*t122*t124;
        double t336 = t3*t6*t11*t61*t87*t124;
        double t147 = -t129;
        double t167 = t15*t121;

        double t185 = -t139;
        double t188 = -t145;

        double t195 = Iwyy*t13*t143;
        double t196 = Iwyy*mt*rf*t143;
        double t197 = Iwyy*t14*t143;
        double t198 = t2*t8*t126;
        double t199 = t7*t8*t126;
        double t200 = t7*t127*dph;
        double t201 = t146/2.0;
        double t202 = -t168;
        double t204 = t190*t190;
        double t209 = t66+t119;
        double t211 = t65+t121;
        double t227 = t5*t189;
        double t228 = t5*t190;
        double t229 = t212*2.0;
        double t230 = t213*2.0;
        double t231 = t213*4.0;
        double t233 = t51+t53+t123;
        double t234 = t57*t173;
        double t235 = t55*t176;
        double t237 = t68*t176;
        double t239 = t107+t146;
        double t240 = Itxx*Iwyy*t11*t186;
        double t241 = Itzz+t212;
        double t242 = t101+t143;
        double t243 = t102+t142;

        double t247 = -t212;
        double t250 = t8*t218;
        double t256 = Itzz*Iwyy*t126*t142;

        double t258 = t246*4.0;
        double t259 = t8*t215;

        double t263 = t68*t177*dth*dps;
        double t267 = rw*t5*t210;
        double t274 = t106+t186;
        double t276 = t5*t246;

        double t279 = t3*t236*4.0;
        double t280 = g*t3*t10*t190*4.0;
        double t282 = rw*t8*t10*t190*4.0;
        double t283 = Iwyy*t87*t126*t145;
        double t286 = t3*t5*t70*t172;
        double t290 = t11*t124*t186;
        double t291 = t57*t61*t172;
        double t292 = t70*t124*t128;

        double t299 = t61*t70*t172*dth;
        double t301 = t288*2.0;
        double t305 = t4*t9*t125*t169;
        double t306 = t288*dph;
        double t307 = t275/2.0;
        double t311 = t16+t149+t153;
        double t312 = t10*t70*t96*t174*dth;
        double t315 = t10*t67*t70*t172;
        double t330 = t6*t11*t100*t124*t126;
        double t331 = t6*t11*t102*t124*t126;
        double t338 = t10*t54*t70*t174*dph;
        double t339 = t51*t287;
        double t346 = t61*t63*t75*t175;
        double t350 = rw*t30*t59*t62*t190;
        double t358 = t215*t236;
        double t368 = t104+t159+t223;
        double t370 = Iwyy+t215+t218;
        double t380 = Itxx+Ityy+t16+t81+t88+t130+t166;
        double t220 = Iwyy*t13*t185;
        double t221 = Iwyy*mt*rf*t185;
        double t222 = Iwyy*t14*t185;

        double t238 = rw*t228;
        double t248 = -t229;
        double t249 = -t231;
        double t251 = t234*2.0;
        double t252 = t235*2.0;
        double t255 = rw*t227*4.0;
        double t264 = -t234;
        double t266 = rw*t5*t209;
        double t268 = rw*t5*t211;
        double t269 = rw*t10*t209;
        double t270 = rw*t211*dth*dom;
        double t271 = Itzz+t247;
        double t272 = t100+t188;
        double t273 = t105+t185;

        double t281 = t235/4.0;
        double t284 = t28+t229;
        double t285 = t33+t230;
        double t289 = -t280;

        double t297 = t54*t239;
        double t298 = t56*t239;
        double t300 = t59*t237*dps;
        double t308 = t276/2.0;

        double t316 = rw*t10*t72*t211;
        double t317 = t3*t10*t241*4.0;
        double t318 = t126+t227;
        double t319 = -t307;
        double t321 = t126+t228;
        double t324 = t147+t203+1.0;
        double t327 = -t299;
        double t329 = rw*t8*t10*t51*t211;
        double t332 = t291/4.0;
        double t334 = -t305;
        double t337 = t53+t123+t202;
        double t343 = -t315;

        double t351 = -t338;
        double t353 = t74*t311;
        double t354 = t190*t242;
        double t356 = -t346;
        double t361 = Iwyy*t190*t243;
        double t366 = rw*t10*t190*t233;

        double t372 = rw*t5*t204*t243;
        double t374 = t228*t290;
        double t375 = rw*t2*t3*t190*t243;
        double t376 = t158+t187+t200;
        double t379 = rw*t7*t8*t190*t243;
        double t387 = t2*t126*t190*t243;
        double t398 = t56*t380;
        double t399 = Itxx+Ityy+t18+t130+t166+t247;
        double t401 = t279+t282;
        double t406 = Ifxx+t80+t87+t137+t138+t239;
        double t409 = t6*t11*t124*t228*t243;
        double t426 = t30*t59*t61*t370;
        double t428 = t51*t124*t244*t274;
        double t438 = t212+t380;
        double t458 = t17+t22+t25+t82+t89+t131+t167+t229;
        double t495 = t160+t161+t162+t195+t196+t197;
        double t663 = t118+t178+t179+t180+t181+t182+t183+t208+t217+t232+t240;
        double t260 = Iwyy+t238;
        double t265 = -t251;
        double t295 = t28+t248;
        double t296 = t34+t249;
        double t303 = t268/4.0;
        double t309 = t285*dps;
        double t314 = t5*t271*dth;

        double t325 = -t297;
        double t326 = -t298;
        double t328 = t67*t268*dph;
        double t335 = t3*t10*t271;
        double t340 = t5*t8*t284;
        double t345 = Iwyy*t8*t10*t87*t238;
        double t359 = g*t8*t318*4.0;
        double t360 = t70*t324;
        double t362 = Iwyy*t2*t3*t321;
        double t363 = t190*t273;
        double t364 = t5*t128*t284;
        double t365 = t3*t184*t268*dph;
        double t369 = t8*t361;
        double t371 = t214+t269;
        double t373 = -t366;
        double t378 = Iwyy+t213+t266;
        double t381 = rw*t59*t71*t321;
        double t382 = rw*t2*t3*t7*t8*t321;
        double t383 = Itzz*t59*t100*t321;
        double t384 = Itzz*t59*t102*t321;
        double t386 = rw*t5*t204*t272;
        double t388 = -t375;
        double t389 = rw*t2*t8*t190*t272;
        double t390 = rw*t3*t7*t190*t272;
        double t393 = t3*t58*t126*t321;
        double t394 = t3*t71*t126*t321;
        double t395 = t70*t124*t337;
        double t396 = t7*t126*t190*t272;
        double t404 = t199+t354;

        double t408 = t399*dom;
        double t410 = t19+t34+t150+t152+t255;
        double t411 = t398/4.0;
        double t412 = t3*t6*t11*t100*t124*t321;
        double t413 = rw*t3*t100*t190*t321;
        double t414 = t3*t6*t11*t102*t124*t321;
        double t415 = t8*t399;
        double t416 = rw*t3*t102*t190*t321;
        double t417 = t53*t401;
        double t419 = t19+t22+t25+t131+t167+t248;
        double t420 = t6*t11*t124*t228*t272;
        double t421 = t59*t101*t190*t321;

        double t423 = t59*t105*t190*t321;
        double t424 = (Iwyy*t399)/2.0;
        double t433 = t69*t406;

        double t453 = t3*t10*t438;
        double t454 = t8*t190*t243*t321;
        double t455 = t69*t438;
        double t461 = (t2*t126*t399)/2.0;
        double t462 = (t7*t126*t399)/2.0;
        double t468 = t3*t7*t190*t243*t321;
        double t470 = t8*t190*t272*t321;
        double t478 = t2*t3*t190*t272*t321;
        double t480 = t267+t399;
        double t482 = t69*t458*dth*dom;
        double t488 = t3*t56*t458*dph*dth;
        double t494 = t3*t56*t184*t438*dph;

        double t509 = (t2*t3*t321*t399)/2.0;
        double t511 = (t3*t7*t321*t399)/2.0;
        double t514 = t3*t56*t458*(t96-dom);
        double t516 = t163+t164+t165+t220+t221+t222;
        double t533 = t236*t495;
        double t563 = Itxx+Ityy+t19+t34+t77+t81+t88+t90+t94+t116+t132+t205+t235+t264;
        double t569 = t190*t272*t495;
        double t606 = (t399*t495)/2.0;
        double t680 = Itxx+Ityy+t16+t31+t35+t130+t148+t166+t235+t268+t291+t398;
        double t691 = t236*t663;
        double t717 = t190*t243*t663;
        double t720 = t2*t3*t321*t663;
        double t721 = t3*t7*t321*t663;
        double t723 = t190*t272*t663;
        double t302 = t8*t260;
        double t323 = t296*dps;
        double t341 = t335*dph;
        double t342 = -t314;
        double t344 = t5*t295*dom;
        double t347 = Iwyy*t72*t260;
        double t352 = t5*t8*t295*2.0;

        double t385 = t5*t11*t103*t124*t260;
        double t391 = t3*t371;

        double t397 = t3*t378*dph*dth;
        double t403 = t246+t335;
        double t405 = -t396;

        double t429 = -t414;
        double t430 = Itzz*t10*t404;
        double t431 = t419*dom;
        double t432 = -t424;

        double t435 = t288+t340;
        double t436 = t415/2.0;

        double t440 = t72*t410;
        double t445 = t3*t419*dph*dth;
        double t448 = t2*t126*t404;
        double t449 = -Itzz*t10*(t198-t363);
        double t450 = Iwyy*t10*t113*t415;
        double t456 = t236+t390;
        double t459 = t59*t433*2.0;

        double t465 = -t461;
        double t466 = -t462;
        double t471 = -t7*t126*(t198-t363);
        double t472 = -t124*dph*(t226-t360);
        double t473 = t236*t404;
        double t474 = t270+t395;
        double t475 = t238*t404;

        double t485 = t214+t433;
        double t489 = t236*(t198-t363);

        double t492 = (t8*t480)/2.0;
        double t493 = t309+t408;
        double t498 = t331+t413;
        double t499 = t3*t7*t321*t404;
        double t500 = t237+t455;

        double t506 = t190*t272*t404;

        double t513 = t455*(t52-t169);
        double t515 = -t190*t243*(t198-t363);

        double t518 = t393+t394;
        double t520 = -t5*(t246-t453);
        double t521 = t306+t327+t364;
        double t522 = Iwyy*t3*t5*t113*t321*t415;

        double t527 = (t7*t126*(t275-t415))/2.0;

        double t529 = rw*t2*t8*t516;
        double t540 = t151+t154+t225+t326+t353;
        double t543 = t236*(t275-t415)*(-1.0/2.0);
        double t546 = (t399*t404)/2.0;
        double t551 = t236*t516;

        double t573 = t3*t7*t321*(t275-t415)*(-1.0/2.0);

        double t577 = t421+t454;
        double t581 = t190*t243*t516;
        double t592 = t22+t25+t78+t79+t82+t89+t91+t92+t95+t134+t135+t219+t252+t265;
        double t600 = (t51*t67*t563)/2.0;

        double t619 = t468+t478;
        double t621 = t409+t509;
        double t622 = ((t198-t363)*(t275-t415))/2.0;
        double t634 = (t399*t516)/2.0;

        double t657 = t238*(t423-t470);

        double t671 = -rw*t3*t7*(t420-t511);
        double t676 = rw*t2*t8*(t420-t511);

        double t694 = ((t275-t415)*(t330-t416))/2.0;

        double t696 = (Itzz*t3*t5*t680)/4.0;
        double t703 = (t7*t126*t680)/4.0;
        double t712 = (Iwyy*t199*t680)/4.0;
        double t713 = (t215*t680)/4.0;

        double t725 = (Iwyy*t8*t238*t680)/4.0;
        double t741 = (t190*t243*t680)/4.0;
        double t744 = (t190*t272*t680)/4.0;
        double t759 = (t399*t680)/8.0;
        double t763 = (t404*t680)/4.0;

        double t792 = t663*(t275-t415)*(-1.0/2.0);

        double t810 = t680*(t275-t415)*(-1.0/8.0);
        double t825 = t109+t110+t111+t112+t113+t114+t115+t171+t201+t281+t303+t332+t381+t411;

        double t400 = t259+t302;
        double t402 = t391/2.0;
        double t418 = t403*dph*dom;
        double t441 = t190*t243*t302;
        double t442 = t435*dph;
        double t443 = t2*t3*t302*t321;
        double t444 = t3*t7*t302*t321;
        double t446 = t190*t272*t302;

        double t464 = -t459;
        double t469 = Iwyy*t456;

        double t486 = t10*t474*2.0;
        double t487 = t302*t404;
        double t491 = t3*t485;

        double t501 = t321*t436;
        double t502 = t3*t493*dph*2.0;

        double t510 = t292+t341+t342;
        double t512 = t8*t500;
        double t526 = t387+t405;
        double t530 = t521*dne*2.0;
        double t531 = Iwyy*t8*t518;
        double t532 = rw*t2*t518;
        double t534 = t383+t430;

        double t539 = rw*t3*t7*t518;
        double t548 = t384+t449;
        double t549 = t290+t492;
        double t552 = t388+t456;
        double t553 = t372+t466;
        double t555 = t236*t518;
        double t557 = t386+t465;
        double t558 = t59*t540*2.0;
        double t559 = t3*t10*t190*t518;
        double t575 = t190*t272*t498;

        double t598 = t448+t471;
        double t599 = t319+t389+t436;
        double t604 = t96*t592;
        double t611 = t391+t520;
        double t620 = t238*t577;
        double t626 = t254+t334+t373+t397+Tomega;
        double t631 = Iwyy*t621;
        double t642 = (t399*t518)/2.0;

        double t652 = rw*t2*t3*t619;
        double t654 = rw*t2*t8*t619;
        double t655 = rw*t3*t7*t619;
        double t656 = rw*t7*t8*t619;

        double t700 = t506+t515;
        double t714 = -Iwyy*(t499+t2*t3*t321*(t198-t363));
        double t749 = t302*(t499+t2*t3*t321*(t198-t363));
        double t802 = t569+t581;
        double t822 = t356+t759;
        double t824 = t412+t741;
        double t827 = Iwyy*t825;
        double t828 = t429+t744;
        double t832 = -t2*t126*(t346-t759);

        double t836 = ((t499+t2*t3*t321*(t198-t363))*(t275-t415))/2.0;
        double t840 = t7*t126*(t346-t759);
        double t842 = -rw*t2*t8*(t414-t744);
        double t867 = t302*(t346-t759);
        double t898 = t551+t721;
        double t930 = t634+t723;

        double t467 = t236*t400;
        double t476 = t250+t400;
        double t497 = t491/2.0;
        double t519 = t510*dne*4.0;

        double t535 = -t530;
        double t536 = rw*t2*t526;
        double t538 = t8*t532;
        double t541 = Iwyy*t534;
        double t542 = dth*(t344-t442)*(-1.0/2.0);
        double t545 = rw*t3*t7*t526;

        double t550 = Iwyy*t548;
        double t556 = Itzz*t10*t549;
        double t560 = rw*t7*t553;
        double t561 = t286+t512;

        double t565 = rw*t2*(t276-t491)*(-1.0/2.0);
        double t566 = t2*t126*t549;
        double t567 = t7*t126*t549;
        double t568 = rw*t2*t557;

        double t571 = t238*t534;
        double t574 = t287+t316+t464;
        double t578 = t3*t5*t87*t552;

        double t583 = t238*t548;
        double t584 = t236*t549;
        double t589 = t8*t321*t526;
        double t590 = t374+t501;
        double t593 = t236*(t276-t491)*(-1.0/2.0);
        double t601 = t302*t549;
        double t603 = Iwyy*t599;
        double t607 = t362+t532;
        double t608 = rw*t2*t598;
        double t616 = rw*t3*t7*t598;
        double t618 = t190*t243*t549;
        double t623 = t2*t3*t321*t549;
        double t624 = t3*t7*t321*t549;
        double t627 = t190*t272*t549;
        double t630 = -t620;
        double t633 = t190*t243*(t276-t491)*(-1.0/2.0);

        double t636 = t3*t7*t321*(t276-t491)*(-1.0/2.0);
        double t637 = t361*(t276-t491)*(-1.0/2.0);
        double t643 = (t2*t3*t321*(t276-t491))/2.0;

        double t649 = t190*t272*(t276-t491)*(-1.0/2.0);
        double t650 = t8*t631;

        double t658 = t236*t598;
        double t662 = t3*t5*t113*t611;
        double t666 = (t2*t126*t611)/2.0;
        double t667 = (t7*t126*t611)/2.0;
        double t668 = (t399*t534)/2.0;
        double t672 = Itzz*t3*t74*t190*t598;
        double t678 = t399*(t276-t491)*(-1.0/4.0);
        double t679 = (t399*t548)/2.0;
        double t681 = (t236*t611)/2.0;
        double t682 = (t238*t611)/2.0;
        double t685 = (t3*t10*t190*t611)/2.0;
        double t689 = -rw*t2*(t446+(t2*t126*(t275-t415))/2.0);

        double t693 = t379+t599;
        double t698 = (t190*t243*t611)/2.0;
        double t708 = (t190*t272*t611)/2.0;

        double t716 = rw*t2*t700;
        double t719 = (t399*t598)/2.0;
        double t728 = ((t275-t415)*(t276-t491))/4.0;

        double t731 = (t399*t611)/4.0;
        double t733 = (t404*t611)/2.0;
        double t735 = -Iwyy*(t382-t402+(t5*(t246-t453))/2.0);

        double t747 = t526*(t276-t491)*(-1.0/2.0);

        double t756 = (t400*t680)/4.0;
        double t760 = (t598*(t275-t415))/2.0;
        double t767 = t323+t431+t604;
        double t813 = (t549*t611)/2.0;

        double t815 = rw*t2*t8*t802;
        double t816 = rw*t3*t7*t802;
        double t817 = rw*t7*t8*t802;
        double t820 = (t611*(t276-t491))/4.0;
        double t830 = (t526*t680)/4.0;
        double t834 = t117+t133+t206+t325+t343+t440+t558;

        double t843 = ((t441+t527)*(t276-t491))/2.0;
        double t849 = t238*t824;
        double t850 = (t549*t680)/4.0;
        double t852 = ((t276-t491)*(t446+(t2*t126*(t275-t415))/2.0))/2.0;
        double t864 = t539+t703;
        double t870 = t531+t714;
        double t896 = (t611*t663)/2.0;
        double t905 = Itzz*t3*t5*(t696-t10*t113*(t276-t491));
        double t916 = rw*t2*t8*t898;
        double t917 = rw*t3*t7*t898;

        double t929 = t680*(t446+(t2*t126*(t275-t415))/2.0)*(-1.0/4.0);
        double t935 = rw*t2*t8*t930;
        double t936 = rw*t3*t7*t930;
        double t971 = -Iwyy*t8*(-t575+t642+t190*t243*(t330-t416));

        double t976 = -Itzz*t3*t5*(-t575+t642+t190*t243*(t330-t416));

        double t979 = -rw*t3*t7*(-t575+t642+t190*t243*(t330-t416));
        double t980 = rw*t2*(-t575+t642+t190*t243*(t330-t416));
        double t982 = t655+t824;

        double t989 = t652+t828;
        double t990 = t671+t822;

        double t1086 = -rw*t2*(t832+t238*(t414-t744));
        double t1090 = Iwyy*t8*(t832+t238*(t414-t744));
        double t1092 = rw*t2*(t832+t238*(t414-t744));

        double t544 = t8*t536;
        double t579 = t52*t561;
        double t582 = t574*dom;

        double t588 = t283+t541;

        double t594 = Itzz*t10*t399*t476*(-1.0/2.0);
        double t595 = t191+t565;
        double t596 = t256+t550;
        double t597 = Iwyy*t590;
        double t609 = t8*t603;
        double t610 = t361+t536;

        double t613 = t2*t126*t590;
        double t614 = t7*t126*t590;
        double t615 = t8*t608;
        double t628 = -t618;
        double t629 = t336+t556;

        double t639 = rw*t2*t3*t607;
        double t640 = rw*t7*t8*t607;
        double t648 = rw*t7*(t331-t545);

        double t665 = t302*t590;
        double t673 = -t667;

        double t701 = Itzz*t10*t693;

        double t710 = -Iwyy*(t308-t497+rw*t3*t7*(t198-t363));
        double t722 = t8*t716;

        double t736 = -rw*t7*(t475-t567);
        double t742 = rw*t2*(t566-t238*(t198-t363));
        double t743 = rw*t7*t8*t87*(t475-t567);
        double t748 = t385+t682;
        double t751 = t473+t623;
        double t754 = t473+t633;
        double t757 = t444+t666;
        double t758 = t3*t74*t87*t190*(t475-t567);
        double t761 = t495+t608;

        double t765 = Itzz*t3*t74*t190*(t566-t238*(t198-t363));

        double t774 = t3*t767;
        double t775 = t489+t649;
        double t778 = t190*t243*(t566-t238*(t198-t363));

        double t797 = t236*(t489-t624);
        double t799 = t3*t10*t190*(t489-t624);
        double t800 = t559+t589;

        double t835 = t253+t333+t418+t428+t542+Tneta;
        double t841 = rw*t7*(t616-(t7*t126*(t276-t491))/2.0);
        double t846 = -t843;
        double t847 = (Iwyy*t834)/4.0;
        double t853 = (Itzz*t10*t834)/4.0;
        double t854 = ((t275-t415)*(t475-t567))/2.0;

        double t859 = (t2*t126*t834)/4.0;
        double t860 = (t7*t126*t834)/4.0;

        double t865 = ((t566-t238*(t198-t363))*(t275-t415))/2.0;
        double t873 = t8*t321*(t627-(t399*(t198-t363))/2.0);
        double t875 = rw*t7*t864;
        double t876 = t584+t678;
        double t879 = (t236*t834)/4.0;
        double t880 = (t238*t834)/4.0;
        double t887 = t432+t560+t568;
        double t888 = rw*t2*t8*t870;
        double t889 = rw*t3*t7*t870;

        double t903 = t238*(t685-(t8*t321*(t275-t415))/2.0);
        double t904 = t573+t708;
        double t906 = ((t275-t415)*(t489-t624))/2.0;
        double t907 = (t190*t243*t834)/4.0;
        double t910 = (t2*t3*t321*t834)/4.0;
        double t911 = (t3*t7*t321*t834)/4.0;
        double t913 = (t190*t272*t834)/4.0;

        double t942 = t658+t747;
        double t947 = t643+t763;
        double t952 = t555+t830;
        double t959 = -t7*t126*(t636+(t680*(t198-t363))/4.0);
        double t985 = t593+t850;
        double t993 = t681+t810;

        double t995 = rw*t7*t8*t982;
        double t996 = Iwyy*t990;
        double t997 = t97+t289+t312+t339+t445+t488+t513+t519;
        double t1000 = t3*t5*t87*t982;
        double t1004 = Itzz*t3*t5*t989;
        double t1011 = t399*(t636+(t680*(t198-t363))/4.0)*(-1.0/2.0);
        double t1015 = (t680*t834)/1.6e+1;
        double t1027 = -Itzz*t10*(-t654+t698+(t2*t3*t321*(t275-t415))/2.0);
        double t1030 = -rw*t2*t3*(-t654+t698+(t2*t3*t321*(t275-t415))/2.0);
        double t1031 = t543+t676+t731;
        double t1033 = ((t636+(t680*(t198-t363))/4.0)*(t275-t415))/2.0;
        double t1073 = t840+t849;
        double t1074 = t631+t980;
        double t1091 = t8*t1086;
        double t1104 = t263+t359+t365+t482+t486+t494+t502+t535+t600;
        double t587 = -t582;
        double t602 = t8*t597;

        double t625 = t2*t126*t595;
        double t638 = Iwyy*t629;
        double t644 = rw*t2*t3*t610;
        double t646 = rw*t7*t8*t610;

        double t661 = t10*t87*t629;
        double t664 = t7*t126*t629;
        double t669 = t2*t126*t629;
        double t697 = t190*t243*t629;
        double t706 = t190*t272*t629;
        double t711 = (t399*t588)/2.0;
        double t718 = (t399*t596)/2.0;
        double t745 = Itzz*t8*t742;
        double t766 = rw*t2*t757;
        double t768 = t2*t126*t751;

        double t771 = rw*t3*t7*t757;
        double t772 = rw*t7*t8*t761;
        double t776 = rw*t2*t3*t761;
        double t777 = Iwyy*t8*(t443+t673);
        double t781 = t236*t751;

        double t783 = t3*t10*t190*t751;
        double t784 = t238*t754;
        double t785 = rw*t2*t8*t775;
        double t788 = t236*t757;
        double t793 = t302*t751;
        double t798 = t190*t243*t748;

        double t804 = t190*t272*t748;
        double t806 = t238*t775;
        double t807 = t546+t628;
        double t809 = t190*t272*t751;
        double t812 = t190*t243*t757;
        double t818 = Itzz*t3*t5*t800;

        double t829 = t441+t527+t544;
        double t845 = (t399*t757)/2.0;

        double t866 = -t860;
        double t871 = t614+t630;
        double t881 = -t880;
        double t884 = t613+t657;

        double t899 = t2*t126*t876;
        double t900 = t7*t126*t876;
        double t902 = t578+t701;
        double t908 = -t903;

        double t912 = t3*t5*t10*t30*t887;
        double t914 = -t910;
        double t920 = rw*t3*t7*t904;
        double t922 = t302*t876;
        double t923 = t443+t538+t673;
        double t938 = t190*t243*(t347-t847);
        double t945 = rw*t2*t942;

        double t949 = Iwyy*t947;
        double t950 = rw*t2*(t859-t302*(t198-t363));
        double t953 = t2*t126*t947;
        double t955 = t328+t351+t514+t774;
        double t958 = t236*(t859-t302*(t198-t363));

        double t961 = rw*t2*t8*t952;
        double t962 = rw*t7*t8*t952;
        double t963 = t238*t947;
        double t965 = t3*t5*t87*t952;
        double t970 = t302*t947;
        double t975 = t190*t243*(t859-t302*(t198-t363));
        double t984 = t190*t272*t947;
        double t987 = Itzz*t985;
        double t988 = Iwyy*t985;
        double t1001 = t8*t996;
        double t1002 = t2*t126*t985;
        double t1003 = t7*t126*t985;
        double t1007 = (t399*(t859-t302*(t198-t363)))/2.0;
        double t1008 = (t399*t947)/2.0;
        double t1013 = t662+t853;
        double t1014 = t238*t993;

        double t1018 = t190*t243*t985;
        double t1019 = t190*t272*t985;
        double t1023 = t656+t904;
        double t1024 = t622+t913;
        double t1036 = Iwyy*t1031;
        double t1042 = t728+t879;
        double t1043 = t985*(t275-t415)*(-1.0/2.0);
        double t1058 = -t7*t126*(t911+(t611*(t198-t363))/2.0);
        double t1064 = -t236*(t911+(t611*(t198-t363))/2.0);
        double t1066 = t799+t873;
        double t1071 = t813+t879;
        double t1076 = t663+t736+t742;
        double t1078 = Iwyy*t8*t1073;
        double t1083 = rw*t2*t3*t1074;
        double t1084 = rw*t7*t8*t1074;
        double t1095 = t399*(t911+(t611*(t198-t363))/2.0)*(-1.0/2.0);
        double t1101 = -rw*t7*(-t722+t907+(t404*(t275-t415))/2.0);
        double t1106 = -rw*t2*(-t719+t778+t190*t272*(t475-t567));
        double t1120 = t820+t1015;
        double t1170 = t1000+t1027;
        double t1173 = t639+t713+t827+t875;
        double t1182 = t979+t1073;

        double t1204 = t10*t87*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567)));
        double t1206 = -rw*t7*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567)));
        double t659 = -t644;
        double t674 = -t664;
        double t675 = t37+t300+t587;
        double t702 = t345+t638;
        double t707 = -t697;
        double t838 = rw*t7*t829;
        double t848 = -t845;
        double t862 = t8*t321*t807;
        double t863 = t583+t669;
        double t882 = Iwyy*t8*t871;
        double t883 = rw*t7*t871;
        double t892 = Itzz*t3*t5*t871;
        double t894 = Iwyy*t8*t884;
        double t895 = rw*t2*t884;
        double t901 = Itzz*t3*t5*t884;

        double t924 = t238*t902;
        double t925 = rw*t7*t923;
        double t933 = t679+t706;
        double t943 = t487+t866;
        double t948 = t8*t945;
        double t957 = t955*dth;
        double t991 = t672+t818;

        double t1006 = -t1003;
        double t1009 = -t1008;
        double t1010 = t601+t881;
        double t1016 = Itzz*t10*t1013;
        double t1020 = t756+t771;
        double t1028 = rw*t2*t1024;
        double t1029 = Itzz*t10*t1023;
        double t1034 = rw*t3*t7*t1024;
        double t1035 = t665+t908;
        double t1040 = t8*t1036;
        double t1044 = t238*t1042;
        double t1046 = t733+t914;
        double t1067 = Iwyy*t1066;
        double t1068 = rw*t2*t1066;

        double t1080 = t2*t126*t1071;
        double t1081 = t7*t126*t1071;
        double t1082 = -rw*t2*(t806-t899);
        double t1093 = t236*t1071;
        double t1094 = t661+t987;
        double t1096 = t190*t243*t1071;
        double t1097 = t190*t272*t1071;
        double t1107 = t8*t1106;

        double t1110 = t852+t958;
        double t1111 = -rw*t2*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624));
        double t1113 = -rw*t3*t7*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624));
        double t1117 = -t3*t10*t190*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624));
        double t1125 = t87*t1120;
        double t1126 = -Itzz*t10*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673));
        double t1127 = -Iwyy*t8*(t798+(t399*(t443+t673))/2.0+(t498*(t275-t415))/2.0);
        double t1130 = -rw*t2*t3*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673));
        double t1134 = t2*t126*t1120;
        double t1135 = t7*t126*t1120;
        double t1141 = t190*t243*t1120;
        double t1142 = t190*t272*t1120;
        double t1146 = ((t275-t415)*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624)))/2.0;
        double t1148 = (t399*t1120)/2.0;
        double t1152 = -rw*t3*t7*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624));
        double t1156 = t953+t959;
        double t1159 = ((t276-t491)*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673)))/2.0;

        double t1163 = Itzz*(t1002-t238*(t636+(t680*(t198-t363))/4.0));
        double t1165 = rw*t2*(t1002-t238*(t636+(t680*(t198-t363))/4.0));
        double t1166 = t625+t710+t776+t841;
        double t1168 = t533+t637+t816+t945;
        double t1180 = rw*t7*t8*(t650+rw*t2*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)));
        double t1183 = rw*t7*t1182;
        double t1196 = -rw*t7*t8*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0));
        double t1211 = -rw*t7*(t798+(t399*(t443+t673))/2.0+t8*t980+(t498*(t275-t415))/2.0);
        double t1213 = -rw*t2*t3*(t798+(t399*(t443+t673))/2.0+t8*t980+(t498*(t275-t415))/2.0);
        double t1219 = t797+t1011+t1019;
        double t1233 = -rw*t2*t3*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)));
        double t1244 = -rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)));
        double t1291 = -Iwyy*(t1120+rw*t2*t8*(t636+(t680*(t198-t363))/4.0)+rw*t3*t7*(t911+(t611*(t198-t363))/2.0));
        double t1300 = -rw*t7*(t961+t236*(t443+t673)+(t680*(t441+t527))/4.0+rw*t3*t7*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673)));
        double t1323 = t842+t920+t993+t995+t1030;

        double t686 = t675*dph*2.0;
        double t737 = t3*t10*t190*t702;
        double t750 = t190*t243*t702;
        double t752 = t190*t272*t702;
        double t844 = -t838;
        double t855 = t571+t674;
        double t872 = rw*t2*t863;
        double t886 = t10*t87*t863;
        double t928 = -t925;
        double t931 = t668+t707;
        double t940 = rw*t2*t933;
        double t954 = t236*t943;
        double t973 = t190*t272*t943;
        double t981 = t358+t469+t648+t659;
        double t999 = (t399*t943)/2.0;
        double t1012 = Itzz*t1010;
        double t1022 = Iwyy*t8*t1020;
        double t1025 = t594+t924;
        double t1032 = t190*t243*t1010;
        double t1038 = t190*t272*t1010;
        double t1045 = t615+t943;
        double t1047 = -t1044;
        double t1048 = Iwyy*t1046;
        double t1052 = t758+t892;
        double t1053 = t2*t126*t1046;

        double t1057 = t765+t901;
        double t1060 = t236*t1046;
        double t1061 = t238*t1046;
        double t1069 = t190*t272*t1046;
        double t1070 = -t1068;
        double t1087 = -t1081;
        double t1088 = (t399*t1046)/2.0;
        double t1089 = t8*t1082;
        double t1099 = t302*t1094;
        double t1112 = t8*t1111;
        double t1114 = rw*t2*t1110;
        double t1116 = Itzz*t1113;
        double t1123 = t694+t804+t848;
        double t1138 = -t1135;
        double t1139 = t597+t883+t895;
        double t1149 = -t1148;

        double t1157 = rw*t2*t1156;
        double t1160 = t963+t1006;
        double t1167 = Itzz*t3*t5*t1166;
        double t1171 = rw*t7*t8*t1168;
        double t1172 = t1004+t1029;
        double t1184 = t533+t720+t1111;
        double t1185 = t8*t1183;
        double t1186 = -t1183;
        double t1208 = t965+t1126;
        double t1210 = t781+t1009+t1018;
        double t1216 = rw*t2*(t1134-t302*(t636+(t680*(t198-t363))/4.0));
        double t1221 = Iwyy*t1219;
        double t1222 = rw*t2*t1219;
        double t1225 = t7*t126*t1219;
        double t1227 = t302*t1219;
        double t1240 = t785+t1034+t1042;
        double t1246 = t8*t1244;
        double t1269 = t190*t243*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624));

        double t1271 = -rw*t7*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567)));
        double t1272 = t905+t1016+t1125;
        double t1275 = t906+t1095+t1097;
        double t1294 = -rw*t3*t7*(Itzz*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624))+t3*t5*t87*t863);
        double t1297 = t788+t929+t962+t1130;
        double t1303 = t1033+t1064+t1142;
        double t1315 = -t7*t126*(-t1093+t1148+(t985*(t275-t415))/2.0);
        double t1316 = t7*t126*(-t1093+t1148+(t985*(t275-t415))/2.0);
        double t1324 = Iwyy*t1323;

        double t868 = Itzz*t10*t855;
        double t869 = rw*t7*t855;
        double t877 = Itzz*t3*t5*t855;
        double t932 = rw*t7*t931;
        double t941 = -t940;
        double t998 = t3*t5*t87*t981;
        double t1039 = -t1032;
        double t1041 = t3*t5*t87*t1025;
        double t1049 = rw*t7*t1045;
        double t1055 = rw*t7*t1052;
        double t1063 = rw*t2*t1057;
        double t1115 = Itzz*t1112;
        double t1128 = rw*t2*t1123;
        double t1129 = Iwyy*t8*t1123;
        double t1137 = rw*t3*t7*t1123;
        double t1140 = t603+t646+t689+t844;
        double t1158 = t8*t1157;

        double t1164 = t87*t1160;
        double t1169 = t350+t743+t745+t1012;
        double t1175 = t640+t735+t766+t928;
        double t1188 = rw*t2*t3*t1184;
        double t1189 = rw*t7*t8*t1184;
        double t1194 = t760+t973+t975;
        double t1212 = Iwyy*t1210;
        double t1215 = t886+t1163;
        double t1217 = t2*t126*t1210;
        double t1218 = -t1216;
        double t1223 = t8*t1222;
        double t1224 = t302*t1210;
        double t1237 = t865+t1007+t1038;
        double t1241 = Iwyy*t1240;
        double t1247 = t749+t1053+t1058;
        double t1253 = t793+t1061+t1087;
        double t1268 = rw*t2*t3*(t836-t1069+t190*t243*(t911+(t611*(t198-t363))/2.0));

        double t1276 = Iwyy*t1275;
        double t1277 = t238*t1272;
        double t1278 = rw*t2*t1275;
        double t1279 = rw*t3*t7*t1275;
        double t1282 = t718+t752+t1206;
        double t1285 = -rw*t7*(t1160+rw*t3*t7*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624)));
        double t1296 = rw*t7*t8*(t712+t889-t949-t1157);
        double t1299 = rw*t2*t1297;
        double t1304 = -t238*(-t1060+t1141+(t947*(t275-t415))/2.0);
        double t1305 = t238*t1303;
        double t1311 = t1043+t1093+t1149;
        double t1321 = t1152+t1210;
        double t1336 = -rw*t7*(-t1088+t1096+(t751*(t275-t415))/2.0+rw*t2*t8*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)));
        double t1337 = t894+t1067+t1271;
        double t1344 = t996+t1083+t1092+t1186;
        double t874 = -t869;
        double t1133 = -t1128;
        double t1144 = Itzz*t10*t1140;
        double t1154 = t8*t321*t1140;
        double t1155 = t450+t932+t941;
        double t1181 = t3*t10*t190*t1175;

        double t1197 = rw*t2*t1194;
        double t1198 = rw*t7*t1194;
        double t1220 = rw*t2*t8*t1215;
        double t1236 = t854+t999+t1039;
        double t1238 = rw*t2*t1237;
        double t1239 = rw*t3*t7*t1237;
        double t1242 = (t680*t1194)/4.0;
        double t1248 = rw*t2*t1247;
        double t1249 = rw*t3*t7*t1247;
        double t1250 = t236*t1247;
        double t1252 = (t399*t1247)/2.0;
        double t1254 = Itzz*t1253;
        double t1258 = t522+t737+t1055+t1063;
        double t1262 = t190*t272*t1253;
        double t1280 = -t1278;
        double t1287 = -Itzz*(t347+t529-t772-t847-t950+t1049);
        double t1306 = t868+t1116+t1164;
        double t1322 = rw*t7*t1321;
        double t1325 = t867+t1014+t1091+t1137;
        double t1327 = t1112+t1253;
        double t1345 = Itzz*t3*t5*t1344;
        double t1352 = -Itzz*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)));
        double t1356 = rw*t7*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)));
        double t1362 = t1196+t1268+t1303;
        double t1364 = t3*t10*t190*(-t725+t917+t988+t1165+t1188+t1285);
        double t1368 = t1223+t1279+t1311;
        double t1382 = t1299+t1300+t1324;
        double t1153 = t702+t872+t874;
        double t1199 = t3*t1198;
        double t1200 = -t1197;

        double t1230 = t998+t1144;
        double t1251 = -t1250;
        double t1260 = t3*t5*t87*t1258;
        double t1263 = -t1262;
        double t1288 = t10*t1287;
        double t1301 = t1154+t1181;
        double t1307 = rw*t7*t8*t1306;
        double t1317 = t1107+t1236;
        double t1326 = Iwyy*t8*t1325;
        double t1329 = rw*t7*t1327;
        double t1338 = t877+t1115+t1254;
        double t1339 = t777+t888+t1048+t1248;
        double t1346 = t1036+t1084+t1133+t1211;
        double t1348 = t970+t1138+t1158+t1249;
        double t1360 = t1185+t1213+t1325;
        double t1367 = t1204+t1352;
        double t1369 = Iwyy*t1368;
        double t1381 = t1040+t1180+t1280+t1336;
        double t1394 = t1090+t1221+t1233+t1356;
        double t1232 = Itzz*t3*t5*t1230;
        double t1313 = rw*t2*t3*(t815+t938+t1200-(t495*(t275-t415))/2.0);
        double t1318 = rw*t7*t1317;
        double t1319 = rw*t2*t3*t1317;

        double t1331 = t846+t948+t954+t1199;
        double t1340 = rw*t2*t3*t1338;
        double t1341 = rw*t2*t3*t1339;
        double t1342 = t1167+t1288;
        double t1347 = -t10*t87*t1346;
        double t1349 = rw*t7*t1348;
        double t1350 = t1159+t1242+t1251;
        double t1372 = t1146+t1252+t1263+t1269;
        double t1384 = -t3*t10*t190*(-t896+t916-t1189+t1329+t236*(t347-t847)+rw*t2*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)));

        double t1332 = rw*t7*t1331;
        double t1374 = rw*t2*t1372;
        double t1375 = rw*t7*t1372;
        double t1386 = t1345+t1347;
        double t1389 = t922+t1047+t1089+t1239+t1246+t1319;
        double t1405 = -Itzz*(t1384+t8*t321*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))));
        double t1406 = t1022+t1218+t1291+t1296+t1341+t1349;
        double t1333 = -t1332;
        double t1376 = t3*t1375;
        double t1390 = Itzz*t1389;
        double t1408 = t1260+t1405;
        double t1396 = Itzz*(t1114+t1171+t1241+t1313+t1333-Iwyy*t8*(t467+rw*t3*t7*(t446+(t2*t126*(t275-t415))/2.0)));
        double t1409 = rw*t8*t1408;
        double t1414 = -rw*t7*(-t1224+t1315+t1376+t238*(-t1060+t1141+(t947*(t275-t415))/2.0)+rw*t2*t8*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0))));
        double t1410 = t1409*4.0;

        double et1 = t1409;
        double et2 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et3 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double t1421 = -1.0/(et1+et2+et3);
        double et4 = t1410;
        double et5 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et6 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double t1422 = -1.0/(et4+et5+et6);

        double et7 = t1410;
        double et8 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et9 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et10 = t1409;
        double et11 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et12 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et13 = t1409;
        double et14 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et15 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et16 = t1409;
        double et17 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et18 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et19 = t1422*(Itzz*t1344+t30*t74*t887+Itzz*rw*t8*t1139)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352)))+t835*t1421*(t1345-Itzz*t10*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+rw*t8*t1258);
        double et20 = t1104*t1422*(t912+Itzz*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))-Itzz*rw*t3*t8*t10*t190*t1076)+(t997*(Itzz*(-t725+t917+t988+t1165+t1188+t1285)+t10*t87*t1153+Itzz*rw*t72*t321*t1076))/(et7+et8+et9);
        double et21 = (t626*(Itzz*(t1001-t1222-t1322+rw*t2*t3*(t650+rw*t2*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))-Itzz*t10*t1155+rw*t8*t87*(-t602+t1068+rw*t7*(t783-t862))))/(et10+et11+et12);
        double et22 = (rw*t368*(-Itzz*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+t10*t87*(-t711+t750+rw*t2*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567))))+Itzz*rw*t8*(t882+Iwyy*(t783-t862)+rw*t2*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))))))/(et13+et14+et15)+rw*t376*t1421*(Itzz*t1394+t10*t87*t1282+Itzz*rw*t8*t1337);
        double et23 = (rw*t3*t52*(Itzz*(t1364+t8*t321*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t10*t1258))/(et16+et17+et18);
        double et24 = t1410;
        double et25 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et26 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et27 = t1409;
        double et28 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et29 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et30 = t1409;
        double et31 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et32 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et33 = t1409;
        double et34 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et35 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et36 = t1409;
        double et37 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et38 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et39 = ((t912-Itzz*t1346)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et24+et25+et26)+(t835*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+t3*t5*t87*t1346))/(et27+et28+et29)+t1104*t1422*(Itzz*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))-t30*t59*t61*t887)+(t626*(Itzz*t1381+Itzz*t3*t5*t1155))/(et30+et31+et32);
        double et40 = t997*t1422*(Itzz*(-t896+t916-t1189+t1329+t236*(t347-t847)+rw*t2*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)))+t3*t5*t87*t1153)+(rw*t368*(Itzz*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))-t3*t5*t87*(-t711+t750+rw*t2*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567))))))/(et33+et34+et35);
        double et41 = (rw*t376*(Itzz*(t1129-t1276-t1375+rw*t7*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+t3*t5*t87*t1282))/(et36+et37+et38)+rw*t3*t52*t1408*t1421;
        double et42 = t1410;
        double et43 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et44 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et45 = t1409;
        double et46 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et47 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et48 = t1409;
        double et49 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et50 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et51 = t87*(-t1093+t1148+(t985*(t275-t415))/2.0)+Itzz*t10*((t399*t1013)/2.0+(t629*(t275-t415))/2.0)-rw*t8*(Itzz*((t590*t834)/4.0-t549*(t685-(t8*t321*(t275-t415))/2.0))+Itzz*t8*t1070-Itzz*t3*t5*(t3*t10*t190*t629+t3*t5*t113*t321*t415)+rw*t7*t8*t87*(t783-t862))-rw*t2*t3*(Itzz*(-t1088+t1096+(t751*(t275-t415))/2.0)-Itzz*t3*t5*t931+Itzz*rw*t2*t8*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))+rw*t3*t7*(Itzz*t1275-t3*t5*t87*t933);
        double et52 = Itzz*t3*t5*(t236*t629+(t399*(t696-t10*t113*(t276-t491)))/2.0)+rw*t7*t8*(Itzz*t1152+Itzz*t1210+Itzz*t10*t931)+rw*t2*t8*(Itzz*t1219+t10*t87*t933);
        double et53 = t997*t1422*(t1099+t1220+t1277+t1294+t1307+t1340+rw*t72*t321*t1169)+((t87*t1360+Itzz*t10*t1025+rw*t8*(Itzz*t1035+Itzz*t8*t883+Itzz*t8*t895-rw*t30*t59*t61*t74*t204))*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et42+et43+et44)+t1104*t1422*(t1041+t1390+rw*t3*t8*t10*t190*t1169)+(t835*(t10*t1390+rw*t8*(t8*t1055+t8*t1063+Itzz*t3*t5*t1035+t3*t74*t190*t1012)+t3*t5*t87*t1360))/(et45+et46+et47)+(t626*(et51+et52))/(et48+et49+et50);
        double et54 = rw*t376*t1421*(t87*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+Itzz*t10*(Itzz*t10*t1237+t3*t5*t87*t1123)+rw*t8*(t87*(t8*t321*t1237+t3*t10*t190*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)))+Itzz*t3*t5*t1057+Itzz*rw*t7*t8*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))))+rw*t7*t8*t1367+rw*t2*t3*(Itzz*t1372-Itzz*t3*t5*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567))))-Itzz*t3*t5*(Itzz*t10*(t806-t899)+t3*t5*t87*(t832+t238*(t414-t744))));
        double et55 = rw*t368*t1421*(Itzz*(t1224+t1304+t1316)-t10*t87*(Itzz*t10*t1236+Itzz*t3*t5*(t798+(t399*(t443+t673))/2.0+(t498*(t275-t415))/2.0))+rw*t8*(-Itzz*(t8*t321*t1236-t3*t10*t190*t1253)+t3*t5*t87*t1052+Itzz*rw*t2*t8*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))))-t3*t5*t87*(Itzz*t10*(t784-t900)-t3*t5*t87*t1073)+rw*t2*t8*t1367-rw*t3*t7*(Itzz*t1372-Itzz*t3*t5*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567)))))+rw*t3*t122*t1421*(t8*t321*(t1041+t1390)-t3*t10*t190*(t1099+t1220+t1277+t1294+t1307+t1340));
        double et56 = t1410;
        double et57 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et58 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et59 = -Itzz*(-t8*t1324+rw*t2*t1362+rw*t7*(-t1060+t1141+(t947*(t275-t415))/2.0+rw*t2*t8*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0))+rw*t3*t7*(t836-t1069+t190*t243*(t911+(t611*(t198-t363))/2.0))));
        double et60 = rw*t8*(Itzz*(t8*t321*(t609-t1028+t1101+rw*t7*t8*(t369-t716))-t3*t10*t190*(Iwyy*t8*(t382-t402+(t5*(t246-t453))/2.0)-rw*t2*(t911+(t611*(t198-t363))/2.0)+rw*t7*(-t733+t910+rw*t2*t8*(t499+t2*t3*t321*(t198-t363)))-rw*t7*t8*(t8*t362+rw*t2*(t499+t2*t3*t321*(t198-t363)))))-t3*t5*t87*(rw*t7*(t3*t10*t190*t534+Itzz*t3*t8*t228*t243*t321)-rw*t2*(t3*t10*t190*t548+t3*t8*t87*t228*t272*t321)+Itzz*Iwyy*t3*t8*t74*t190));
        double et61 = t10*t87*(Itzz*t10*(t609-t1028+t1101+rw*t7*t8*(t369-t716))-Itzz*t3*t5*(t8*t469+rw*t7*(t754+rw*t3*t7*t700)+rw*t2*t775-rw*t2*t3*(t369-t716)))+Itzz*t3*t5*(rw*t2*t1172-rw*t7*t1170);
        double et62 = t1409;
        double et63 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et64 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et65 = t1410;
        double et66 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et67 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et68 = Itzz*(Iwyy*(-t1060+t1141+(t947*(t275-t415))/2.0+rw*t2*t8*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0))+rw*t3*t7*(t836-t1069+t190*t243*(t911+(t611*(t198-t363))/2.0)))+rw*t2*t1350-Iwyy*t8*(t961+t236*(t443+t673)+(t680*(t441+t527))/4.0+rw*t3*t7*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673))))-rw*t8*(Itzz*(t8*t321*(t815+t938+t1200-(t495*(t275-t415))/2.0)+t3*t10*t190*t1339)+Itzz*t3*t5*(rw*t2*t991+t3*t10*t190*t588+Itzz*Iwyy*t3*t8*t228*t243*t321));
        double et69 = Itzz*t10*(Itzz*t10*(t815+t938+t1200-(t495*(t275-t415))/2.0)+t3*t5*t87*t1168)+Itzz*t3*t5*(Iwyy*t1170+rw*t2*t1208);
        double et70 = t1409;
        double et71 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et72 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et73 = rw*t376*(Itzz*(Iwyy*t1362+Iwyy*t8*t1297-rw*t7*t1350)-t10*t87*(Itzz*t10*(-t817+t1198+(t516*(t275-t415))/2.0+t190*t272*(t347-t847))-t3*t5*t87*(t551+rw*t7*t942-rw*t2*t3*t802+(Iwyy*t190*t272*(t276-t491))/2.0))-rw*t8*(Itzz*(t8*t321*(-t817+t1198+(t516*(t275-t415))/2.0+t190*t272*(t347-t847))-t3*t10*t190*(-Iwyy*(t911+(t611*(t198-t363))/2.0)+Iwyy*t8*t757+rw*t7*t1247+rw*t7*t8*t870))-Itzz*t3*t5*(rw*t7*t991+t3*t10*t190*t596+Iwyy*t3*t8*t87*t228*t272*t321))+t3*t5*t87*(Iwyy*t1172+rw*t7*t1208));
        double et74 = 1.0/(et70+et71+et72);
        double et75 = (t997*(Itzz*t1406-t10*t87*t1342+rw*t72*t321*(t426+Itzz*(t347+t529-t772-t847-t950+t1049))+Itzz*t3*t5*(Itzz*t10*t1175-t3*t5*t87*t1173)))/(et56+et57+et58)+t626*t1421*(et59+et60+et61)+t1104*t1422*(t1232+t1396-rw*t3*t8*t10*t190*(t426+Itzz*(t347+t529-t772-t847-t950+t1049)))+(t835*(t10*t1396+rw*t8*(Itzz*t3*t5*t1301+t3*t74*t190*t1287)+Itzz*t3*t5*t1382))/(et62+et63+et64)+((Itzz*t1382+rw*t8*(Itzz*t1301+t30*t59*t74*t228*t370)+t10*t87*t1230)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et65+et66+et67)+rw*t368*t1421*(et68+et69);
        double et76 = et73*et74+rw*t3*t122*t1421*(t8*t321*(t1232+t1396)+t3*t10*t190*(Itzz*t1406-t10*t87*t1342+Itzz*t3*t5*(Itzz*t10*t1175-t3*t5*t87*t1173)));
        double et77 = t1410;
        double et78 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et79 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);

        double et80 = t1410;
        double et81 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et82 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et83 = t1409;
        double et84 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et85 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et86 = t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))-rw*t8*(t1384+t8*t321*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))));
        double et87 = rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))));
        double et88 = t1409;
        double et89 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et90 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et91 = t10*t87*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+Itzz*t3*t5*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))));
        double et92 = -Itzz*rw*t3*t5*t8*(t882+Iwyy*(t783-t862)+rw*t2*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))));
        double et93 = t1409;
        double et94 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et95 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et96 = t1409;
        double et97 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et98 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et99 = t1409;
        double et100 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et101 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et102 = ((t1386+Itzz*rw*t3*t5*t8*t1139)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et77+et78+et79);
        double et103 = (t1104*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))-Itzz*t8*t10*t59*t238*t1076))/(et80+et81+et82)+(t835*(et86+et87))/(et83+et84+et85);
        double et104 = (t626*(-Itzz*t10*t1381+t3*t5*t87*(t1001-t1222-t1322+rw*t2*t3*(t650+rw*t2*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+Itzz*rw*t3*t5*t8*(-t602+t1068+rw*t7*(t783-t862))))/(et88+et89+et90)+t997*t1422*(t10*t87*(-t896+t916-t1189+t1329+t236*(t347-t847)+rw*t2*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)))+Itzz*t3*t5*(-t725+t917+t988+t1165+t1188+t1285)+Itzz*rw*t3*t5*t72*t321*t1076)+(rw*t368*(et91+et92))/(et93+et94+et95);
        double et105 = (rw*t376*(t10*t87*(t1129-t1276-t1375+rw*t7*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+Itzz*t3*t5*t1394+Itzz*rw*t3*t5*t8*t1337))/(et96+et97+et98);
        double et106 = (rw*t3*t122*(Itzz*t10*(t1384+t8*t321*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))+Itzz*t3*t5*(t1364+t8*t321*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))))/(et99+et100+et101);

        d_dOM[index] = et75+et76;
        
        index = index + num_threads;
        // index1 = (grid_size[0] > 1) ? index : 0;
        index2 = (grid_size[1] > 1) ? index : 0;
        index3 = (grid_size[2] > 1) ? index : 0;
        index4 = (grid_size[3] > 1) ? index : 0;
        index5 = (grid_size[4] > 1) ? index : 0;
        index6 = (grid_size[5] > 1) ? index : 0;
        index7 = (grid_size[6] > 1) ? index : 0;
        index8 = (grid_size[7] > 1) ? index : 0;

        uindex1 = (active_actions[0] == 1) ? index : 0;
        uindex2 = (active_actions[1] == 1) ? index : 0;
    }
}

__global__ void dyn7_mex_continuous(double* const d_dNE, // outputs
    double* const PH, double* const TH, double* const OM, double* const dPH, double* const dTH,
    double* const dPS, double* const dOM, double* const dNE, // input states
    double const* const in1, double const* const in2, // input actions
    const double mw, const double Iwxx, const double Iwyy, const double Iwzz, 
    const double mf, const double Ifxx, const double Ifyy, const double Ifzz, 
    const double mt, const double Itxx, const double Ityy, const double Itzz, 
    const double nw, const double nt, const double rw, const double rf, const double rt, const double fcoeff, const double g, 
    int32_t const* const grid_size, int32_t const* const active_actions)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3]
                        * grid_size[4] * grid_size[5] * grid_size[6] * grid_size[7];
    double th, ps, om, ne, dph, dth, dps, dom, dne, Tneta, Tomega;
    ps = 0; ne = 0;

    // int index1 = (grid_size[0] > 1) ? index : 0;
    int index2 = (grid_size[1] > 1) ? index : 0;
    int index3 = (grid_size[2] > 1) ? index : 0;
    int index4 = (grid_size[3] > 1) ? index : 0;
    int index5 = (grid_size[4] > 1) ? index : 0;
    int index6 = (grid_size[5] > 1) ? index : 0;
    int index7 = (grid_size[6] > 1) ? index : 0;
    int index8 = (grid_size[7] > 1) ? index : 0;

    int uindex1 = (active_actions[0] == 1) ? index : 0;
    int uindex2 = (active_actions[1] == 1) ? index : 0;

    while (index < num_elements)
    {
        th = TH[index2];
        om = OM[index3];
        dph = dPH[index4];
        dth = dTH[index5];
        dps = dPS[index6];
        dom = dOM[index7];
        dne = dNE[index8];

        Tomega = nw*in1[uindex1];
        Tneta = nt*in2[uindex2];

        double t2 = cos(ph);
        double t3 = cos(th);
        double t4 = cos(ps);
        double t5 = cos(om);
        double t6 = cos(ne);
        double t7 = sin(ph);
        double t8 = sin(th);
        double t9 = sin(ps);
        double t10 = sin(om);
        double t11 = sin(ne);
        double t12 = mf+mt;
        double t13 = mf*rf;
        double t14 = mt*rt;
        double t15 = rf+rt;
        double t16 = Ifxx*2.0;
        double t17 = Ifxx*4.0;
        double t18 = Ifyy*2.0;
        double t19 = Ifyy*4.0;
        double t20 = Ifzz*2.0;
        double t21 = Ifzz*4.0;
        double t22 = Itxx*2.0;
        double t23 = Itxx*3.0;
        double t24 = Itxx*4.0;
        double t25 = Ityy*2.0;
        double t26 = Ityy*3.0;
        double t27 = Ityy*4.0;
        double t28 = Itzz*2.0;
        double t29 = Itzz*4.0;
        double t30 = Itzz*Itzz;
        double t31 = Iwxx*2.0;
        double t32 = Iwxx*4.0;
        double t33 = Iwyy*2.0;
        double t34 = Iwyy*4.0;
        double t35 = Iwzz*2.0;
        double t36 = Iwzz*4.0;
        double t37 = fcoeff*2.0;

        double t40 = mt*2.0;
        double t42 = rf*rf;

        double t44 = rw*rw;
        double t45 = Tomega*4.0;
        double t46 = th*2.0;
        double t47 = ps*2.0;
        double t48 = om*2.0;
        double t49 = ne*2.0;
        double t50 = dom*2.0;
        double t51 = dph*dph;
        double t52 = dth*dth;
        double t53 = dom*dom;
        double t79 = Ifyy*8.0;
        double t80 = -Ifzz;
        double t83 = -Ityy;
        double t87 = -Itzz;
        double t92 = Iwyy*8.0;
        double t93 = -Iwzz;
        double t109 = Ifxx/2.0;
        double t110 = Ifzz/2.0;
        double t111 = Itxx/4.0;
        double t112 = Ityy/4.0;
        double t113 = Itzz/2.0;
        double t114 = Iwxx/2.0;
        double t115 = Iwzz/2.0;
        double t54 = cos(t46);
        double t55 = cos(t47);
        double t56 = cos(t48);
        double t57 = cos(t49);
        double t58 = t2*t2;
        double t59 = t3*t3;
        double t60 = t4*t4;
        double t61 = t5*t5;
        double t62 = t5*t5*t5;
        double t63 = t6*t6;
        double t64 = t13*2.0;
        double t65 = t13*4.0;
        double t66 = t14*2.0;
        double t67 = sin(t46);
        double t68 = sin(t47);
        double t69 = sin(t48);
        double t70 = sin(t49);
        double t71 = t7*t7;
        double t72 = t8*t8;
        double t73 = t9*t9;
        double t74 = t10*t10;
        double t75 = t11*t11;
        double t76 = mw+t12;
        double t77 = -t16;
        double t78 = -t17;
        double t81 = -t20;
        double t82 = -t21;
        double t84 = -t25;
        double t85 = -t26;
        double t86 = -t27;
        double t88 = -t28;
        double t89 = -t29;
        double t90 = -t31;
        double t91 = -t32;
        double t94 = -t35;
        double t95 = -t36;
        double t96 = t8*dph;
        double t97 = -t45;
        double t98 = rf*t12;
        double t99 = mt*t15;
        double t100 = t2*t5;
        double t101 = t2*t10;
        double t102 = t5*t7;
        double t103 = t6*t8;
        double t104 = t2*dph*dps;
        double t105 = t7*t10;
        double t106 = t8*t11;
        double t107 = rf*t13;
        double t108 = t15*t15;
        double t116 = rf*t14*4.0;
        double t117 = rf*t14*6.0;
        double t118 = Ifyy*Iwyy*t8;
        double t120 = t15*t40;
        double t122 = -t52;
        double t124 = Itxx+t83;
        double t125 = Iwxx+t93;
        double t133 = rt*t14*3.0;
        double t134 = rt*t14*4.0;
        double t135 = rf*t14*8.0;
        double t140 = t3*t6*t10;
        double t141 = t51+t52;
        double t144 = t3*t10*t11;
        double t148 = t20+t28;
        double t155 = t12*2.0;
        double t156 = t12*3.0;
        double t157 = t2*t3*dph*dth*2.0;
        double t158 = t7*t8*t52;
        double t159 = t3*t7*dph*dth*2.0;
        double t179 = Iwyy*mt*t8*t42;
        double t180 = Iwyy*rt*t8*t14;
        double t181 = (Itxx*Iwyy*t8)/2.0;
        double t182 = (Ityy*Iwyy*t8)/2.0;
        double t183 = rf*t8*t14*t33;
        double t191 = Iwyy*rw*t3*t7*t8;
        double t219 = t12*t42*4.0;
        double t119 = t98*2.0;
        double t121 = t99*4.0;
        double t123 = t50*t96;
        double t126 = rw*t76;
        double t127 = t96+dps;
        double t128 = t96+dom;
        double t129 = t54*3.0;
        double t130 = rf*t64;
        double t131 = rf*t65;
        double t132 = rt*t66;

        double t137 = Itxx*t63;
        double t138 = Ityy*t75;
        double t139 = t8*t100;
        double t142 = t8*t101;
        double t143 = t8*t102;
        double t145 = t8*t105;
        double t146 = t15*t99;
        double t149 = t22*t63;
        double t150 = t27*t63;
        double t151 = t35*t60;
        double t152 = t24*t75;
        double t153 = t25*t75;
        double t154 = t31*t73;
        double t160 = Iwyy*t13*t101;
        double t161 = Iwyy*mt*rf*t101;
        double t162 = Iwyy*t14*t101;
        double t163 = Iwyy*t13*t105;
        double t164 = Iwyy*mt*rf*t105;
        double t165 = Iwyy*t14*t105;
        double t166 = t40*t108;
        double t168 = t51*t54;
        double t169 = t51*t59;
        double t171 = t107/2.0;
        double t172 = t22+t84;
        double t173 = t23+t85;
        double t174 = t24+t86;
        double t175 = t124*t124;
        double t176 = t31+t94;
        double t177 = t32+t95;
        double t178 = Iwyy*t8*t107;
        double t184 = t50+t96;
        double t186 = -t140;
        double t187 = -t157;
        double t189 = t14+t98;
        double t190 = t13+t99;
        double t203 = t56*t59*2.0;
        double t205 = t42*t155;
        double t206 = t42*t156;

        double t208 = Ityy*Iwyy*t11*t140;
        double t210 = t64+t120;
        double t212 = t57*t124;
        double t213 = t55*t125;
        double t214 = t68*t125;
        double t215 = t44*t58*t76;

        double t217 = t57*t182;
        double t218 = t44*t71*t76;
        double t223 = t2*t8*t141;
        double t225 = t61*t148;
        double t226 = t10*t57*t67*4.0;
        double t232 = Itxx*Iwyy*t8*t57*(-1.0/2.0);
        double t236 = t5*t6*t11*t124;
        double t244 = t103+t144;
        double t246 = t8*t70*t124;
        double t253 = t6*t11*t53*t124;
        double t254 = t4*t9*t52*t125;

        double t275 = t3*t10*t70*t124;
        double t287 = t5*t67*t70*t124;
        double t288 = t3*t69*t70*t124;

        double t333 = t6*t11*t61*t122*t124;
        double t336 = t3*t6*t11*t61*t87*t124;
        double t147 = -t129;
        double t167 = t15*t121;

        double t185 = -t139;
        double t188 = -t145;

        double t195 = Iwyy*t13*t143;
        double t196 = Iwyy*mt*rf*t143;
        double t197 = Iwyy*t14*t143;
        double t198 = t2*t8*t126;
        double t199 = t7*t8*t126;
        double t200 = t7*t127*dph;
        double t201 = t146/2.0;
        double t202 = -t168;
        double t204 = t190*t190;
        double t209 = t66+t119;
        double t211 = t65+t121;
        double t227 = t5*t189;
        double t228 = t5*t190;
        double t229 = t212*2.0;
        double t230 = t213*2.0;
        double t231 = t213*4.0;
        double t233 = t51+t53+t123;
        double t234 = t57*t173;
        double t235 = t55*t176;
        double t237 = t68*t176;
        double t239 = t107+t146;
        double t240 = Itxx*Iwyy*t11*t186;
        double t241 = Itzz+t212;
        double t242 = t101+t143;
        double t243 = t102+t142;

        double t247 = -t212;
        double t250 = t8*t218;
        double t256 = Itzz*Iwyy*t126*t142;

        double t258 = t246*4.0;
        double t259 = t8*t215;

        double t263 = t68*t177*dth*dps;
        double t267 = rw*t5*t210;
        double t274 = t106+t186;
        double t276 = t5*t246;

        double t279 = t3*t236*4.0;
        double t280 = g*t3*t10*t190*4.0;
        double t282 = rw*t8*t10*t190*4.0;
        double t283 = Iwyy*t87*t126*t145;
        double t286 = t3*t5*t70*t172;
        double t290 = t11*t124*t186;
        double t291 = t57*t61*t172;
        double t292 = t70*t124*t128;

        double t299 = t61*t70*t172*dth;
        double t301 = t288*2.0;
        double t305 = t4*t9*t125*t169;
        double t306 = t288*dph;
        double t307 = t275/2.0;
        double t311 = t16+t149+t153;
        double t312 = t10*t70*t96*t174*dth;
        double t315 = t10*t67*t70*t172;
        double t330 = t6*t11*t100*t124*t126;
        double t331 = t6*t11*t102*t124*t126;
        double t338 = t10*t54*t70*t174*dph;
        double t339 = t51*t287;
        double t346 = t61*t63*t75*t175;
        double t350 = rw*t30*t59*t62*t190;
        double t358 = t215*t236;
        double t368 = t104+t159+t223;
        double t370 = Iwyy+t215+t218;
        double t380 = Itxx+Ityy+t16+t81+t88+t130+t166;
        double t220 = Iwyy*t13*t185;
        double t221 = Iwyy*mt*rf*t185;
        double t222 = Iwyy*t14*t185;

        double t238 = rw*t228;
        double t248 = -t229;
        double t249 = -t231;
        double t251 = t234*2.0;
        double t252 = t235*2.0;
        double t255 = rw*t227*4.0;
        double t264 = -t234;
        double t266 = rw*t5*t209;
        double t268 = rw*t5*t211;
        double t269 = rw*t10*t209;
        double t270 = rw*t211*dth*dom;
        double t271 = Itzz+t247;
        double t272 = t100+t188;
        double t273 = t105+t185;

        double t281 = t235/4.0;
        double t284 = t28+t229;
        double t285 = t33+t230;
        double t289 = -t280;

        double t297 = t54*t239;
        double t298 = t56*t239;
        double t300 = t59*t237*dps;
        double t308 = t276/2.0;

        double t316 = rw*t10*t72*t211;
        double t317 = t3*t10*t241*4.0;
        double t318 = t126+t227;
        double t319 = -t307;
        double t321 = t126+t228;
        double t324 = t147+t203+1.0;
        double t327 = -t299;
        double t329 = rw*t8*t10*t51*t211;
        double t332 = t291/4.0;
        double t334 = -t305;
        double t337 = t53+t123+t202;
        double t343 = -t315;

        double t351 = -t338;
        double t353 = t74*t311;
        double t354 = t190*t242;
        double t356 = -t346;
        double t361 = Iwyy*t190*t243;
        double t366 = rw*t10*t190*t233;

        double t372 = rw*t5*t204*t243;
        double t374 = t228*t290;
        double t375 = rw*t2*t3*t190*t243;
        double t376 = t158+t187+t200;
        double t379 = rw*t7*t8*t190*t243;
        double t387 = t2*t126*t190*t243;
        double t398 = t56*t380;
        double t399 = Itxx+Ityy+t18+t130+t166+t247;
        double t401 = t279+t282;
        double t406 = Ifxx+t80+t87+t137+t138+t239;
        double t409 = t6*t11*t124*t228*t243;
        double t426 = t30*t59*t61*t370;
        double t428 = t51*t124*t244*t274;
        double t438 = t212+t380;
        double t458 = t17+t22+t25+t82+t89+t131+t167+t229;
        double t495 = t160+t161+t162+t195+t196+t197;
        double t663 = t118+t178+t179+t180+t181+t182+t183+t208+t217+t232+t240;
        double t260 = Iwyy+t238;
        double t265 = -t251;
        double t295 = t28+t248;
        double t296 = t34+t249;
        double t303 = t268/4.0;
        double t309 = t285*dps;
        double t314 = t5*t271*dth;

        double t325 = -t297;
        double t326 = -t298;
        double t328 = t67*t268*dph;
        double t335 = t3*t10*t271;
        double t340 = t5*t8*t284;
        double t345 = Iwyy*t8*t10*t87*t238;
        double t359 = g*t8*t318*4.0;
        double t360 = t70*t324;
        double t362 = Iwyy*t2*t3*t321;
        double t363 = t190*t273;
        double t364 = t5*t128*t284;
        double t365 = t3*t184*t268*dph;
        double t369 = t8*t361;
        double t371 = t214+t269;
        double t373 = -t366;
        double t378 = Iwyy+t213+t266;
        double t381 = rw*t59*t71*t321;
        double t382 = rw*t2*t3*t7*t8*t321;
        double t383 = Itzz*t59*t100*t321;
        double t384 = Itzz*t59*t102*t321;
        double t386 = rw*t5*t204*t272;
        double t388 = -t375;
        double t389 = rw*t2*t8*t190*t272;
        double t390 = rw*t3*t7*t190*t272;
        double t393 = t3*t58*t126*t321;
        double t394 = t3*t71*t126*t321;
        double t395 = t70*t124*t337;
        double t396 = t7*t126*t190*t272;
        double t404 = t199+t354;

        double t408 = t399*dom;
        double t410 = t19+t34+t150+t152+t255;
        double t411 = t398/4.0;
        double t412 = t3*t6*t11*t100*t124*t321;
        double t413 = rw*t3*t100*t190*t321;
        double t414 = t3*t6*t11*t102*t124*t321;
        double t415 = t8*t399;
        double t416 = rw*t3*t102*t190*t321;
        double t417 = t53*t401;
        double t419 = t19+t22+t25+t131+t167+t248;
        double t420 = t6*t11*t124*t228*t272;
        double t421 = t59*t101*t190*t321;

        double t423 = t59*t105*t190*t321;
        double t424 = (Iwyy*t399)/2.0;
        double t433 = t69*t406;

        double t453 = t3*t10*t438;
        double t454 = t8*t190*t243*t321;
        double t455 = t69*t438;
        double t461 = (t2*t126*t399)/2.0;
        double t462 = (t7*t126*t399)/2.0;
        double t468 = t3*t7*t190*t243*t321;
        double t470 = t8*t190*t272*t321;
        double t478 = t2*t3*t190*t272*t321;
        double t480 = t267+t399;
        double t482 = t69*t458*dth*dom;
        double t488 = t3*t56*t458*dph*dth;
        double t494 = t3*t56*t184*t438*dph;

        double t509 = (t2*t3*t321*t399)/2.0;
        double t511 = (t3*t7*t321*t399)/2.0;
        double t514 = t3*t56*t458*(t96-dom);
        double t516 = t163+t164+t165+t220+t221+t222;
        double t533 = t236*t495;
        double t563 = Itxx+Ityy+t19+t34+t77+t81+t88+t90+t94+t116+t132+t205+t235+t264;
        double t569 = t190*t272*t495;
        double t606 = (t399*t495)/2.0;
        double t680 = Itxx+Ityy+t16+t31+t35+t130+t148+t166+t235+t268+t291+t398;
        double t691 = t236*t663;
        double t717 = t190*t243*t663;
        double t720 = t2*t3*t321*t663;
        double t721 = t3*t7*t321*t663;
        double t723 = t190*t272*t663;
        double t302 = t8*t260;
        double t323 = t296*dps;
        double t341 = t335*dph;
        double t342 = -t314;
        double t344 = t5*t295*dom;
        double t347 = Iwyy*t72*t260;
        double t352 = t5*t8*t295*2.0;

        double t385 = t5*t11*t103*t124*t260;
        double t391 = t3*t371;

        double t397 = t3*t378*dph*dth;
        double t403 = t246+t335;
        double t405 = -t396;

        double t429 = -t414;
        double t430 = Itzz*t10*t404;
        double t431 = t419*dom;
        double t432 = -t424;

        double t435 = t288+t340;
        double t436 = t415/2.0;

        double t440 = t72*t410;
        double t445 = t3*t419*dph*dth;
        double t448 = t2*t126*t404;
        double t449 = -Itzz*t10*(t198-t363);
        double t450 = Iwyy*t10*t113*t415;
        double t456 = t236+t390;
        double t459 = t59*t433*2.0;

        double t465 = -t461;
        double t466 = -t462;
        double t471 = -t7*t126*(t198-t363);
        double t472 = -t124*dph*(t226-t360);
        double t473 = t236*t404;
        double t474 = t270+t395;
        double t475 = t238*t404;

        double t485 = t214+t433;
        double t489 = t236*(t198-t363);

        double t492 = (t8*t480)/2.0;
        double t493 = t309+t408;
        double t498 = t331+t413;
        double t499 = t3*t7*t321*t404;
        double t500 = t237+t455;

        double t506 = t190*t272*t404;

        double t513 = t455*(t52-t169);
        double t515 = -t190*t243*(t198-t363);

        double t518 = t393+t394;
        double t520 = -t5*(t246-t453);
        double t521 = t306+t327+t364;
        double t522 = Iwyy*t3*t5*t113*t321*t415;

        double t527 = (t7*t126*(t275-t415))/2.0;

        double t529 = rw*t2*t8*t516;
        double t540 = t151+t154+t225+t326+t353;
        double t543 = t236*(t275-t415)*(-1.0/2.0);
        double t546 = (t399*t404)/2.0;
        double t551 = t236*t516;

        double t573 = t3*t7*t321*(t275-t415)*(-1.0/2.0);

        double t577 = t421+t454;
        double t581 = t190*t243*t516;
        double t592 = t22+t25+t78+t79+t82+t89+t91+t92+t95+t134+t135+t219+t252+t265;
        double t600 = (t51*t67*t563)/2.0;

        double t619 = t468+t478;
        double t621 = t409+t509;
        double t622 = ((t198-t363)*(t275-t415))/2.0;
        double t634 = (t399*t516)/2.0;

        double t657 = t238*(t423-t470);

        double t671 = -rw*t3*t7*(t420-t511);
        double t676 = rw*t2*t8*(t420-t511);

        double t694 = ((t275-t415)*(t330-t416))/2.0;

        double t696 = (Itzz*t3*t5*t680)/4.0;
        double t703 = (t7*t126*t680)/4.0;
        double t712 = (Iwyy*t199*t680)/4.0;
        double t713 = (t215*t680)/4.0;

        double t725 = (Iwyy*t8*t238*t680)/4.0;
        double t741 = (t190*t243*t680)/4.0;
        double t744 = (t190*t272*t680)/4.0;
        double t759 = (t399*t680)/8.0;
        double t763 = (t404*t680)/4.0;

        double t792 = t663*(t275-t415)*(-1.0/2.0);

        double t810 = t680*(t275-t415)*(-1.0/8.0);
        double t825 = t109+t110+t111+t112+t113+t114+t115+t171+t201+t281+t303+t332+t381+t411;

        double t400 = t259+t302;
        double t402 = t391/2.0;
        double t418 = t403*dph*dom;
        double t441 = t190*t243*t302;
        double t442 = t435*dph;
        double t443 = t2*t3*t302*t321;
        double t444 = t3*t7*t302*t321;
        double t446 = t190*t272*t302;

        double t464 = -t459;
        double t469 = Iwyy*t456;

        double t486 = t10*t474*2.0;
        double t487 = t302*t404;
        double t491 = t3*t485;

        double t501 = t321*t436;
        double t502 = t3*t493*dph*2.0;

        double t510 = t292+t341+t342;
        double t512 = t8*t500;
        double t526 = t387+t405;
        double t530 = t521*dne*2.0;
        double t531 = Iwyy*t8*t518;
        double t532 = rw*t2*t518;
        double t534 = t383+t430;

        double t539 = rw*t3*t7*t518;
        double t548 = t384+t449;
        double t549 = t290+t492;
        double t552 = t388+t456;
        double t553 = t372+t466;
        double t555 = t236*t518;
        double t557 = t386+t465;
        double t558 = t59*t540*2.0;
        double t559 = t3*t10*t190*t518;
        double t575 = t190*t272*t498;

        double t598 = t448+t471;
        double t599 = t319+t389+t436;
        double t604 = t96*t592;
        double t611 = t391+t520;
        double t620 = t238*t577;
        double t626 = t254+t334+t373+t397+Tomega;
        double t631 = Iwyy*t621;
        double t642 = (t399*t518)/2.0;

        double t652 = rw*t2*t3*t619;
        double t654 = rw*t2*t8*t619;
        double t655 = rw*t3*t7*t619;
        double t656 = rw*t7*t8*t619;

        double t700 = t506+t515;
        double t714 = -Iwyy*(t499+t2*t3*t321*(t198-t363));
        double t749 = t302*(t499+t2*t3*t321*(t198-t363));
        double t802 = t569+t581;
        double t822 = t356+t759;
        double t824 = t412+t741;
        double t827 = Iwyy*t825;
        double t828 = t429+t744;
        double t832 = -t2*t126*(t346-t759);

        double t836 = ((t499+t2*t3*t321*(t198-t363))*(t275-t415))/2.0;
        double t840 = t7*t126*(t346-t759);
        double t842 = -rw*t2*t8*(t414-t744);
        double t867 = t302*(t346-t759);
        double t898 = t551+t721;
        double t930 = t634+t723;

        double t467 = t236*t400;
        double t476 = t250+t400;
        double t497 = t491/2.0;
        double t519 = t510*dne*4.0;

        double t535 = -t530;
        double t536 = rw*t2*t526;
        double t538 = t8*t532;
        double t541 = Iwyy*t534;
        double t542 = dth*(t344-t442)*(-1.0/2.0);
        double t545 = rw*t3*t7*t526;

        double t550 = Iwyy*t548;
        double t556 = Itzz*t10*t549;
        double t560 = rw*t7*t553;
        double t561 = t286+t512;

        double t565 = rw*t2*(t276-t491)*(-1.0/2.0);
        double t566 = t2*t126*t549;
        double t567 = t7*t126*t549;
        double t568 = rw*t2*t557;

        double t571 = t238*t534;
        double t574 = t287+t316+t464;
        double t578 = t3*t5*t87*t552;

        double t583 = t238*t548;
        double t584 = t236*t549;
        double t589 = t8*t321*t526;
        double t590 = t374+t501;
        double t593 = t236*(t276-t491)*(-1.0/2.0);
        double t601 = t302*t549;
        double t603 = Iwyy*t599;
        double t607 = t362+t532;
        double t608 = rw*t2*t598;
        double t616 = rw*t3*t7*t598;
        double t618 = t190*t243*t549;
        double t623 = t2*t3*t321*t549;
        double t624 = t3*t7*t321*t549;
        double t627 = t190*t272*t549;
        double t630 = -t620;
        double t633 = t190*t243*(t276-t491)*(-1.0/2.0);

        double t636 = t3*t7*t321*(t276-t491)*(-1.0/2.0);
        double t637 = t361*(t276-t491)*(-1.0/2.0);
        double t643 = (t2*t3*t321*(t276-t491))/2.0;

        double t649 = t190*t272*(t276-t491)*(-1.0/2.0);
        double t650 = t8*t631;

        double t658 = t236*t598;
        double t662 = t3*t5*t113*t611;
        double t666 = (t2*t126*t611)/2.0;
        double t667 = (t7*t126*t611)/2.0;
        double t668 = (t399*t534)/2.0;
        double t672 = Itzz*t3*t74*t190*t598;
        double t678 = t399*(t276-t491)*(-1.0/4.0);
        double t679 = (t399*t548)/2.0;
        double t681 = (t236*t611)/2.0;
        double t682 = (t238*t611)/2.0;
        double t685 = (t3*t10*t190*t611)/2.0;
        double t689 = -rw*t2*(t446+(t2*t126*(t275-t415))/2.0);

        double t693 = t379+t599;
        double t698 = (t190*t243*t611)/2.0;
        double t708 = (t190*t272*t611)/2.0;

        double t716 = rw*t2*t700;
        double t719 = (t399*t598)/2.0;
        double t728 = ((t275-t415)*(t276-t491))/4.0;

        double t731 = (t399*t611)/4.0;
        double t733 = (t404*t611)/2.0;
        double t735 = -Iwyy*(t382-t402+(t5*(t246-t453))/2.0);

        double t747 = t526*(t276-t491)*(-1.0/2.0);

        double t756 = (t400*t680)/4.0;
        double t760 = (t598*(t275-t415))/2.0;
        double t767 = t323+t431+t604;
        double t813 = (t549*t611)/2.0;

        double t815 = rw*t2*t8*t802;
        double t816 = rw*t3*t7*t802;
        double t817 = rw*t7*t8*t802;
        double t820 = (t611*(t276-t491))/4.0;
        double t830 = (t526*t680)/4.0;
        double t834 = t117+t133+t206+t325+t343+t440+t558;

        double t843 = ((t441+t527)*(t276-t491))/2.0;
        double t849 = t238*t824;
        double t850 = (t549*t680)/4.0;
        double t852 = ((t276-t491)*(t446+(t2*t126*(t275-t415))/2.0))/2.0;
        double t864 = t539+t703;
        double t870 = t531+t714;
        double t896 = (t611*t663)/2.0;
        double t905 = Itzz*t3*t5*(t696-t10*t113*(t276-t491));
        double t916 = rw*t2*t8*t898;
        double t917 = rw*t3*t7*t898;

        double t929 = t680*(t446+(t2*t126*(t275-t415))/2.0)*(-1.0/4.0);
        double t935 = rw*t2*t8*t930;
        double t936 = rw*t3*t7*t930;
        double t971 = -Iwyy*t8*(-t575+t642+t190*t243*(t330-t416));

        double t976 = -Itzz*t3*t5*(-t575+t642+t190*t243*(t330-t416));

        double t979 = -rw*t3*t7*(-t575+t642+t190*t243*(t330-t416));
        double t980 = rw*t2*(-t575+t642+t190*t243*(t330-t416));
        double t982 = t655+t824;

        double t989 = t652+t828;
        double t990 = t671+t822;

        double t1086 = -rw*t2*(t832+t238*(t414-t744));
        double t1090 = Iwyy*t8*(t832+t238*(t414-t744));
        double t1092 = rw*t2*(t832+t238*(t414-t744));

        double t544 = t8*t536;
        double t579 = t52*t561;
        double t582 = t574*dom;

        double t588 = t283+t541;

        double t594 = Itzz*t10*t399*t476*(-1.0/2.0);
        double t595 = t191+t565;
        double t596 = t256+t550;
        double t597 = Iwyy*t590;
        double t609 = t8*t603;
        double t610 = t361+t536;

        double t613 = t2*t126*t590;
        double t614 = t7*t126*t590;
        double t615 = t8*t608;
        double t628 = -t618;
        double t629 = t336+t556;

        double t639 = rw*t2*t3*t607;
        double t640 = rw*t7*t8*t607;
        double t648 = rw*t7*(t331-t545);

        double t665 = t302*t590;
        double t673 = -t667;

        double t701 = Itzz*t10*t693;

        double t710 = -Iwyy*(t308-t497+rw*t3*t7*(t198-t363));
        double t722 = t8*t716;

        double t736 = -rw*t7*(t475-t567);
        double t742 = rw*t2*(t566-t238*(t198-t363));
        double t743 = rw*t7*t8*t87*(t475-t567);
        double t748 = t385+t682;
        double t751 = t473+t623;
        double t754 = t473+t633;
        double t757 = t444+t666;
        double t758 = t3*t74*t87*t190*(t475-t567);
        double t761 = t495+t608;

        double t765 = Itzz*t3*t74*t190*(t566-t238*(t198-t363));

        double t774 = t3*t767;
        double t775 = t489+t649;
        double t778 = t190*t243*(t566-t238*(t198-t363));

        double t797 = t236*(t489-t624);
        double t799 = t3*t10*t190*(t489-t624);
        double t800 = t559+t589;

        double t835 = t253+t333+t418+t428+t542+Tneta;
        double t841 = rw*t7*(t616-(t7*t126*(t276-t491))/2.0);
        double t846 = -t843;
        double t847 = (Iwyy*t834)/4.0;
        double t853 = (Itzz*t10*t834)/4.0;
        double t854 = ((t275-t415)*(t475-t567))/2.0;

        double t859 = (t2*t126*t834)/4.0;
        double t860 = (t7*t126*t834)/4.0;

        double t865 = ((t566-t238*(t198-t363))*(t275-t415))/2.0;
        double t873 = t8*t321*(t627-(t399*(t198-t363))/2.0);
        double t875 = rw*t7*t864;
        double t876 = t584+t678;
        double t879 = (t236*t834)/4.0;
        double t880 = (t238*t834)/4.0;
        double t887 = t432+t560+t568;
        double t888 = rw*t2*t8*t870;
        double t889 = rw*t3*t7*t870;

        double t903 = t238*(t685-(t8*t321*(t275-t415))/2.0);
        double t904 = t573+t708;
        double t906 = ((t275-t415)*(t489-t624))/2.0;
        double t907 = (t190*t243*t834)/4.0;
        double t910 = (t2*t3*t321*t834)/4.0;
        double t911 = (t3*t7*t321*t834)/4.0;
        double t913 = (t190*t272*t834)/4.0;

        double t942 = t658+t747;
        double t947 = t643+t763;
        double t952 = t555+t830;
        double t959 = -t7*t126*(t636+(t680*(t198-t363))/4.0);
        double t985 = t593+t850;
        double t993 = t681+t810;

        double t995 = rw*t7*t8*t982;
        double t996 = Iwyy*t990;
        double t997 = t97+t289+t312+t339+t445+t488+t513+t519;
        double t1000 = t3*t5*t87*t982;
        double t1004 = Itzz*t3*t5*t989;
        double t1011 = t399*(t636+(t680*(t198-t363))/4.0)*(-1.0/2.0);
        double t1015 = (t680*t834)/1.6e+1;
        double t1027 = -Itzz*t10*(-t654+t698+(t2*t3*t321*(t275-t415))/2.0);
        double t1030 = -rw*t2*t3*(-t654+t698+(t2*t3*t321*(t275-t415))/2.0);
        double t1031 = t543+t676+t731;
        double t1033 = ((t636+(t680*(t198-t363))/4.0)*(t275-t415))/2.0;
        double t1073 = t840+t849;
        double t1074 = t631+t980;
        double t1091 = t8*t1086;
        double t1104 = t263+t359+t365+t482+t486+t494+t502+t535+t600;
        double t587 = -t582;
        double t602 = t8*t597;

        double t625 = t2*t126*t595;
        double t638 = Iwyy*t629;
        double t644 = rw*t2*t3*t610;
        double t646 = rw*t7*t8*t610;

        double t661 = t10*t87*t629;
        double t664 = t7*t126*t629;
        double t669 = t2*t126*t629;
        double t697 = t190*t243*t629;
        double t706 = t190*t272*t629;
        double t711 = (t399*t588)/2.0;
        double t718 = (t399*t596)/2.0;
        double t745 = Itzz*t8*t742;
        double t766 = rw*t2*t757;
        double t768 = t2*t126*t751;

        double t771 = rw*t3*t7*t757;
        double t772 = rw*t7*t8*t761;
        double t776 = rw*t2*t3*t761;
        double t777 = Iwyy*t8*(t443+t673);
        double t781 = t236*t751;

        double t783 = t3*t10*t190*t751;
        double t784 = t238*t754;
        double t785 = rw*t2*t8*t775;
        double t788 = t236*t757;
        double t793 = t302*t751;
        double t798 = t190*t243*t748;

        double t804 = t190*t272*t748;
        double t806 = t238*t775;
        double t807 = t546+t628;
        double t809 = t190*t272*t751;
        double t812 = t190*t243*t757;
        double t818 = Itzz*t3*t5*t800;

        double t829 = t441+t527+t544;
        double t845 = (t399*t757)/2.0;

        double t866 = -t860;
        double t871 = t614+t630;
        double t881 = -t880;
        double t884 = t613+t657;

        double t899 = t2*t126*t876;
        double t900 = t7*t126*t876;
        double t902 = t578+t701;
        double t908 = -t903;

        double t912 = t3*t5*t10*t30*t887;
        double t914 = -t910;
        double t920 = rw*t3*t7*t904;
        double t922 = t302*t876;
        double t923 = t443+t538+t673;
        double t938 = t190*t243*(t347-t847);
        double t945 = rw*t2*t942;

        double t949 = Iwyy*t947;
        double t950 = rw*t2*(t859-t302*(t198-t363));
        double t953 = t2*t126*t947;
        double t955 = t328+t351+t514+t774;
        double t958 = t236*(t859-t302*(t198-t363));

        double t961 = rw*t2*t8*t952;
        double t962 = rw*t7*t8*t952;
        double t963 = t238*t947;
        double t965 = t3*t5*t87*t952;
        double t970 = t302*t947;
        double t975 = t190*t243*(t859-t302*(t198-t363));
        double t984 = t190*t272*t947;
        double t987 = Itzz*t985;
        double t988 = Iwyy*t985;
        double t1001 = t8*t996;
        double t1002 = t2*t126*t985;
        double t1003 = t7*t126*t985;
        double t1007 = (t399*(t859-t302*(t198-t363)))/2.0;
        double t1008 = (t399*t947)/2.0;
        double t1013 = t662+t853;
        double t1014 = t238*t993;

        double t1018 = t190*t243*t985;
        double t1019 = t190*t272*t985;
        double t1023 = t656+t904;
        double t1024 = t622+t913;
        double t1036 = Iwyy*t1031;
        double t1042 = t728+t879;
        double t1043 = t985*(t275-t415)*(-1.0/2.0);
        double t1058 = -t7*t126*(t911+(t611*(t198-t363))/2.0);
        double t1064 = -t236*(t911+(t611*(t198-t363))/2.0);
        double t1066 = t799+t873;
        double t1071 = t813+t879;
        double t1076 = t663+t736+t742;
        double t1078 = Iwyy*t8*t1073;
        double t1083 = rw*t2*t3*t1074;
        double t1084 = rw*t7*t8*t1074;
        double t1095 = t399*(t911+(t611*(t198-t363))/2.0)*(-1.0/2.0);
        double t1101 = -rw*t7*(-t722+t907+(t404*(t275-t415))/2.0);
        double t1106 = -rw*t2*(-t719+t778+t190*t272*(t475-t567));
        double t1120 = t820+t1015;
        double t1170 = t1000+t1027;
        double t1173 = t639+t713+t827+t875;
        double t1182 = t979+t1073;

        double t1204 = t10*t87*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567)));
        double t1206 = -rw*t7*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567)));
        double t659 = -t644;
        double t674 = -t664;
        double t675 = t37+t300+t587;
        double t702 = t345+t638;
        double t707 = -t697;
        double t838 = rw*t7*t829;
        double t848 = -t845;
        double t862 = t8*t321*t807;
        double t863 = t583+t669;
        double t882 = Iwyy*t8*t871;
        double t883 = rw*t7*t871;
        double t892 = Itzz*t3*t5*t871;
        double t894 = Iwyy*t8*t884;
        double t895 = rw*t2*t884;
        double t901 = Itzz*t3*t5*t884;

        double t924 = t238*t902;
        double t925 = rw*t7*t923;
        double t933 = t679+t706;
        double t943 = t487+t866;
        double t948 = t8*t945;
        double t957 = t955*dth;
        double t991 = t672+t818;

        double t1006 = -t1003;
        double t1009 = -t1008;
        double t1010 = t601+t881;
        double t1016 = Itzz*t10*t1013;
        double t1020 = t756+t771;
        double t1028 = rw*t2*t1024;
        double t1029 = Itzz*t10*t1023;
        double t1034 = rw*t3*t7*t1024;
        double t1035 = t665+t908;
        double t1040 = t8*t1036;
        double t1044 = t238*t1042;
        double t1046 = t733+t914;
        double t1067 = Iwyy*t1066;
        double t1068 = rw*t2*t1066;

        double t1080 = t2*t126*t1071;
        double t1081 = t7*t126*t1071;
        double t1082 = -rw*t2*(t806-t899);
        double t1093 = t236*t1071;
        double t1094 = t661+t987;
        double t1096 = t190*t243*t1071;
        double t1097 = t190*t272*t1071;
        double t1107 = t8*t1106;

        double t1110 = t852+t958;
        double t1111 = -rw*t2*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624));
        double t1113 = -rw*t3*t7*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624));
        double t1117 = -t3*t10*t190*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624));
        double t1125 = t87*t1120;
        double t1126 = -Itzz*t10*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673));
        double t1127 = -Iwyy*t8*(t798+(t399*(t443+t673))/2.0+(t498*(t275-t415))/2.0);
        double t1130 = -rw*t2*t3*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673));
        double t1134 = t2*t126*t1120;
        double t1135 = t7*t126*t1120;
        double t1141 = t190*t243*t1120;
        double t1142 = t190*t272*t1120;
        double t1146 = ((t275-t415)*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624)))/2.0;
        double t1148 = (t399*t1120)/2.0;
        double t1152 = -rw*t3*t7*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624));
        double t1156 = t953+t959;
        double t1159 = ((t276-t491)*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673)))/2.0;

        double t1163 = Itzz*(t1002-t238*(t636+(t680*(t198-t363))/4.0));
        double t1165 = rw*t2*(t1002-t238*(t636+(t680*(t198-t363))/4.0));
        double t1166 = t625+t710+t776+t841;
        double t1168 = t533+t637+t816+t945;
        double t1180 = rw*t7*t8*(t650+rw*t2*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)));
        double t1183 = rw*t7*t1182;
        double t1196 = -rw*t7*t8*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0));
        double t1211 = -rw*t7*(t798+(t399*(t443+t673))/2.0+t8*t980+(t498*(t275-t415))/2.0);
        double t1213 = -rw*t2*t3*(t798+(t399*(t443+t673))/2.0+t8*t980+(t498*(t275-t415))/2.0);
        double t1219 = t797+t1011+t1019;
        double t1233 = -rw*t2*t3*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)));
        double t1244 = -rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)));
        double t1291 = -Iwyy*(t1120+rw*t2*t8*(t636+(t680*(t198-t363))/4.0)+rw*t3*t7*(t911+(t611*(t198-t363))/2.0));
        double t1300 = -rw*t7*(t961+t236*(t443+t673)+(t680*(t441+t527))/4.0+rw*t3*t7*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673)));
        double t1323 = t842+t920+t993+t995+t1030;

        double t686 = t675*dph*2.0;
        double t737 = t3*t10*t190*t702;
        double t750 = t190*t243*t702;
        double t752 = t190*t272*t702;
        double t844 = -t838;
        double t855 = t571+t674;
        double t872 = rw*t2*t863;
        double t886 = t10*t87*t863;
        double t928 = -t925;
        double t931 = t668+t707;
        double t940 = rw*t2*t933;
        double t954 = t236*t943;
        double t973 = t190*t272*t943;
        double t981 = t358+t469+t648+t659;
        double t999 = (t399*t943)/2.0;
        double t1012 = Itzz*t1010;
        double t1022 = Iwyy*t8*t1020;
        double t1025 = t594+t924;
        double t1032 = t190*t243*t1010;
        double t1038 = t190*t272*t1010;
        double t1045 = t615+t943;
        double t1047 = -t1044;
        double t1048 = Iwyy*t1046;
        double t1052 = t758+t892;
        double t1053 = t2*t126*t1046;

        double t1057 = t765+t901;
        double t1060 = t236*t1046;
        double t1061 = t238*t1046;
        double t1069 = t190*t272*t1046;
        double t1070 = -t1068;
        double t1087 = -t1081;
        double t1088 = (t399*t1046)/2.0;
        double t1089 = t8*t1082;
        double t1099 = t302*t1094;
        double t1112 = t8*t1111;
        double t1114 = rw*t2*t1110;
        double t1116 = Itzz*t1113;
        double t1123 = t694+t804+t848;
        double t1138 = -t1135;
        double t1139 = t597+t883+t895;
        double t1149 = -t1148;

        double t1157 = rw*t2*t1156;
        double t1160 = t963+t1006;
        double t1167 = Itzz*t3*t5*t1166;
        double t1171 = rw*t7*t8*t1168;
        double t1172 = t1004+t1029;
        double t1184 = t533+t720+t1111;
        double t1185 = t8*t1183;
        double t1186 = -t1183;
        double t1208 = t965+t1126;
        double t1210 = t781+t1009+t1018;
        double t1216 = rw*t2*(t1134-t302*(t636+(t680*(t198-t363))/4.0));
        double t1221 = Iwyy*t1219;
        double t1222 = rw*t2*t1219;
        double t1225 = t7*t126*t1219;
        double t1227 = t302*t1219;
        double t1240 = t785+t1034+t1042;
        double t1246 = t8*t1244;
        double t1269 = t190*t243*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624));

        double t1271 = -rw*t7*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567)));
        double t1272 = t905+t1016+t1125;
        double t1275 = t906+t1095+t1097;
        double t1294 = -rw*t3*t7*(Itzz*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624))+t3*t5*t87*t863);
        double t1297 = t788+t929+t962+t1130;
        double t1303 = t1033+t1064+t1142;
        double t1315 = -t7*t126*(-t1093+t1148+(t985*(t275-t415))/2.0);
        double t1316 = t7*t126*(-t1093+t1148+(t985*(t275-t415))/2.0);
        double t1324 = Iwyy*t1323;

        double t868 = Itzz*t10*t855;
        double t869 = rw*t7*t855;
        double t877 = Itzz*t3*t5*t855;
        double t932 = rw*t7*t931;
        double t941 = -t940;
        double t998 = t3*t5*t87*t981;
        double t1039 = -t1032;
        double t1041 = t3*t5*t87*t1025;
        double t1049 = rw*t7*t1045;
        double t1055 = rw*t7*t1052;
        double t1063 = rw*t2*t1057;
        double t1115 = Itzz*t1112;
        double t1128 = rw*t2*t1123;
        double t1129 = Iwyy*t8*t1123;
        double t1137 = rw*t3*t7*t1123;
        double t1140 = t603+t646+t689+t844;
        double t1158 = t8*t1157;

        double t1164 = t87*t1160;
        double t1169 = t350+t743+t745+t1012;
        double t1175 = t640+t735+t766+t928;
        double t1188 = rw*t2*t3*t1184;
        double t1189 = rw*t7*t8*t1184;
        double t1194 = t760+t973+t975;
        double t1212 = Iwyy*t1210;
        double t1215 = t886+t1163;
        double t1217 = t2*t126*t1210;
        double t1218 = -t1216;
        double t1223 = t8*t1222;
        double t1224 = t302*t1210;
        double t1237 = t865+t1007+t1038;
        double t1241 = Iwyy*t1240;
        double t1247 = t749+t1053+t1058;
        double t1253 = t793+t1061+t1087;
        double t1268 = rw*t2*t3*(t836-t1069+t190*t243*(t911+(t611*(t198-t363))/2.0));

        double t1276 = Iwyy*t1275;
        double t1277 = t238*t1272;
        double t1278 = rw*t2*t1275;
        double t1279 = rw*t3*t7*t1275;
        double t1282 = t718+t752+t1206;
        double t1285 = -rw*t7*(t1160+rw*t3*t7*(-t768+t238*(t499+t2*t3*t321*(t198-t363))+t7*t126*(t489-t624)));
        double t1296 = rw*t7*t8*(t712+t889-t949-t1157);
        double t1299 = rw*t2*t1297;
        double t1304 = -t238*(-t1060+t1141+(t947*(t275-t415))/2.0);
        double t1305 = t238*t1303;
        double t1311 = t1043+t1093+t1149;
        double t1321 = t1152+t1210;
        double t1336 = -rw*t7*(-t1088+t1096+(t751*(t275-t415))/2.0+rw*t2*t8*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)));
        double t1337 = t894+t1067+t1271;
        double t1344 = t996+t1083+t1092+t1186;
        double t874 = -t869;
        double t1133 = -t1128;
        double t1144 = Itzz*t10*t1140;
        double t1154 = t8*t321*t1140;
        double t1155 = t450+t932+t941;
        double t1181 = t3*t10*t190*t1175;

        double t1197 = rw*t2*t1194;
        double t1198 = rw*t7*t1194;
        double t1220 = rw*t2*t8*t1215;
        double t1236 = t854+t999+t1039;
        double t1238 = rw*t2*t1237;
        double t1239 = rw*t3*t7*t1237;
        double t1242 = (t680*t1194)/4.0;
        double t1248 = rw*t2*t1247;
        double t1249 = rw*t3*t7*t1247;
        double t1250 = t236*t1247;
        double t1252 = (t399*t1247)/2.0;
        double t1254 = Itzz*t1253;
        double t1258 = t522+t737+t1055+t1063;
        double t1262 = t190*t272*t1253;
        double t1280 = -t1278;
        double t1287 = -Itzz*(t347+t529-t772-t847-t950+t1049);
        double t1306 = t868+t1116+t1164;
        double t1322 = rw*t7*t1321;
        double t1325 = t867+t1014+t1091+t1137;
        double t1327 = t1112+t1253;
        double t1345 = Itzz*t3*t5*t1344;
        double t1352 = -Itzz*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)));
        double t1356 = rw*t7*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)));
        double t1362 = t1196+t1268+t1303;
        double t1364 = t3*t10*t190*(-t725+t917+t988+t1165+t1188+t1285);
        double t1368 = t1223+t1279+t1311;
        double t1382 = t1299+t1300+t1324;
        double t1153 = t702+t872+t874;
        double t1199 = t3*t1198;
        double t1200 = -t1197;

        double t1230 = t998+t1144;
        double t1251 = -t1250;
        double t1260 = t3*t5*t87*t1258;
        double t1263 = -t1262;
        double t1288 = t10*t1287;
        double t1301 = t1154+t1181;
        double t1307 = rw*t7*t8*t1306;
        double t1317 = t1107+t1236;
        double t1326 = Iwyy*t8*t1325;
        double t1329 = rw*t7*t1327;
        double t1338 = t877+t1115+t1254;
        double t1339 = t777+t888+t1048+t1248;
        double t1346 = t1036+t1084+t1133+t1211;
        double t1348 = t970+t1138+t1158+t1249;
        double t1360 = t1185+t1213+t1325;
        double t1367 = t1204+t1352;
        double t1369 = Iwyy*t1368;
        double t1381 = t1040+t1180+t1280+t1336;
        double t1394 = t1090+t1221+t1233+t1356;
        double t1232 = Itzz*t3*t5*t1230;
        double t1313 = rw*t2*t3*(t815+t938+t1200-(t495*(t275-t415))/2.0);
        double t1318 = rw*t7*t1317;
        double t1319 = rw*t2*t3*t1317;

        double t1331 = t846+t948+t954+t1199;
        double t1340 = rw*t2*t3*t1338;
        double t1341 = rw*t2*t3*t1339;
        double t1342 = t1167+t1288;
        double t1347 = -t10*t87*t1346;
        double t1349 = rw*t7*t1348;
        double t1350 = t1159+t1242+t1251;
        double t1372 = t1146+t1252+t1263+t1269;
        double t1384 = -t3*t10*t190*(-t896+t916-t1189+t1329+t236*(t347-t847)+rw*t2*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)));

        double t1332 = rw*t7*t1331;
        double t1374 = rw*t2*t1372;
        double t1375 = rw*t7*t1372;
        double t1386 = t1345+t1347;
        double t1389 = t922+t1047+t1089+t1239+t1246+t1319;
        double t1405 = -Itzz*(t1384+t8*t321*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))));
        double t1406 = t1022+t1218+t1291+t1296+t1341+t1349;
        double t1333 = -t1332;
        double t1376 = t3*t1375;
        double t1390 = Itzz*t1389;
        double t1408 = t1260+t1405;
        double t1396 = Itzz*(t1114+t1171+t1241+t1313+t1333-Iwyy*t8*(t467+rw*t3*t7*(t446+(t2*t126*(t275-t415))/2.0)));
        double t1409 = rw*t8*t1408;
        double t1414 = -rw*t7*(-t1224+t1315+t1376+t238*(-t1060+t1141+(t947*(t275-t415))/2.0)+rw*t2*t8*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0))));
        double t1410 = t1409*4.0;

        double et1 = t1409;
        double et2 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et3 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double t1421 = -1.0/(et1+et2+et3);
        double et4 = t1410;
        double et5 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et6 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double t1422 = -1.0/(et4+et5+et6);

        double et7 = t1410;
        double et8 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et9 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et10 = t1409;
        double et11 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et12 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et13 = t1409;
        double et14 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et15 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et16 = t1409;
        double et17 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et18 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et19 = t1422*(Itzz*t1344+t30*t74*t887+Itzz*rw*t8*t1139)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352)))+t835*t1421*(t1345-Itzz*t10*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+rw*t8*t1258);
        double et20 = t1104*t1422*(t912+Itzz*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))-Itzz*rw*t3*t8*t10*t190*t1076)+(t997*(Itzz*(-t725+t917+t988+t1165+t1188+t1285)+t10*t87*t1153+Itzz*rw*t72*t321*t1076))/(et7+et8+et9);
        double et21 = (t626*(Itzz*(t1001-t1222-t1322+rw*t2*t3*(t650+rw*t2*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))-Itzz*t10*t1155+rw*t8*t87*(-t602+t1068+rw*t7*(t783-t862))))/(et10+et11+et12);
        double et22 = (rw*t368*(-Itzz*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+t10*t87*(-t711+t750+rw*t2*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567))))+Itzz*rw*t8*(t882+Iwyy*(t783-t862)+rw*t2*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))))))/(et13+et14+et15)+rw*t376*t1421*(Itzz*t1394+t10*t87*t1282+Itzz*rw*t8*t1337);
        double et23 = (rw*t3*t52*(Itzz*(t1364+t8*t321*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t10*t1258))/(et16+et17+et18);
        double et24 = t1410;
        double et25 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et26 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et27 = t1409;
        double et28 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et29 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et30 = t1409;
        double et31 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et32 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et33 = t1409;
        double et34 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et35 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et36 = t1409;
        double et37 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et38 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et39 = ((t912-Itzz*t1346)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et24+et25+et26)+(t835*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+t3*t5*t87*t1346))/(et27+et28+et29)+t1104*t1422*(Itzz*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))-t30*t59*t61*t887)+(t626*(Itzz*t1381+Itzz*t3*t5*t1155))/(et30+et31+et32);
        double et40 = t997*t1422*(Itzz*(-t896+t916-t1189+t1329+t236*(t347-t847)+rw*t2*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)))+t3*t5*t87*t1153)+(rw*t368*(Itzz*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))-t3*t5*t87*(-t711+t750+rw*t2*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567))))))/(et33+et34+et35);
        double et41 = (rw*t376*(Itzz*(t1129-t1276-t1375+rw*t7*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+t3*t5*t87*t1282))/(et36+et37+et38)+rw*t3*t52*t1408*t1421;
        double et42 = t1410;
        double et43 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et44 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et45 = t1409;
        double et46 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et47 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et48 = t1409;
        double et49 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et50 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et51 = t87*(-t1093+t1148+(t985*(t275-t415))/2.0)+Itzz*t10*((t399*t1013)/2.0+(t629*(t275-t415))/2.0)-rw*t8*(Itzz*((t590*t834)/4.0-t549*(t685-(t8*t321*(t275-t415))/2.0))+Itzz*t8*t1070-Itzz*t3*t5*(t3*t10*t190*t629+t3*t5*t113*t321*t415)+rw*t7*t8*t87*(t783-t862))-rw*t2*t3*(Itzz*(-t1088+t1096+(t751*(t275-t415))/2.0)-Itzz*t3*t5*t931+Itzz*rw*t2*t8*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))+rw*t3*t7*(Itzz*t1275-t3*t5*t87*t933);
        double et52 = Itzz*t3*t5*(t236*t629+(t399*(t696-t10*t113*(t276-t491)))/2.0)+rw*t7*t8*(Itzz*t1152+Itzz*t1210+Itzz*t10*t931)+rw*t2*t8*(Itzz*t1219+t10*t87*t933);
        double et53 = t997*t1422*(t1099+t1220+t1277+t1294+t1307+t1340+rw*t72*t321*t1169)+((t87*t1360+Itzz*t10*t1025+rw*t8*(Itzz*t1035+Itzz*t8*t883+Itzz*t8*t895-rw*t30*t59*t61*t74*t204))*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et42+et43+et44)+t1104*t1422*(t1041+t1390+rw*t3*t8*t10*t190*t1169)+(t835*(t10*t1390+rw*t8*(t8*t1055+t8*t1063+Itzz*t3*t5*t1035+t3*t74*t190*t1012)+t3*t5*t87*t1360))/(et45+et46+et47)+(t626*(et51+et52))/(et48+et49+et50);
        double et54 = rw*t376*t1421*(t87*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+Itzz*t10*(Itzz*t10*t1237+t3*t5*t87*t1123)+rw*t8*(t87*(t8*t321*t1237+t3*t10*t190*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)))+Itzz*t3*t5*t1057+Itzz*rw*t7*t8*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))))+rw*t7*t8*t1367+rw*t2*t3*(Itzz*t1372-Itzz*t3*t5*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567))))-Itzz*t3*t5*(Itzz*t10*(t806-t899)+t3*t5*t87*(t832+t238*(t414-t744))));
        double et55 = rw*t368*t1421*(Itzz*(t1224+t1304+t1316)-t10*t87*(Itzz*t10*t1236+Itzz*t3*t5*(t798+(t399*(t443+t673))/2.0+(t498*(t275-t415))/2.0))+rw*t8*(-Itzz*(t8*t321*t1236-t3*t10*t190*t1253)+t3*t5*t87*t1052+Itzz*rw*t2*t8*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))))-t3*t5*t87*(Itzz*t10*(t784-t900)-t3*t5*t87*t1073)+rw*t2*t8*t1367-rw*t3*t7*(Itzz*t1372-Itzz*t3*t5*(t976+Itzz*t10*(-t719+t778+t190*t272*(t475-t567)))))+rw*t3*t122*t1421*(t8*t321*(t1041+t1390)-t3*t10*t190*(t1099+t1220+t1277+t1294+t1307+t1340));
        double et56 = t1410;
        double et57 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et58 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et59 = -Itzz*(-t8*t1324+rw*t2*t1362+rw*t7*(-t1060+t1141+(t947*(t275-t415))/2.0+rw*t2*t8*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0))+rw*t3*t7*(t836-t1069+t190*t243*(t911+(t611*(t198-t363))/2.0))));
        double et60 = rw*t8*(Itzz*(t8*t321*(t609-t1028+t1101+rw*t7*t8*(t369-t716))-t3*t10*t190*(Iwyy*t8*(t382-t402+(t5*(t246-t453))/2.0)-rw*t2*(t911+(t611*(t198-t363))/2.0)+rw*t7*(-t733+t910+rw*t2*t8*(t499+t2*t3*t321*(t198-t363)))-rw*t7*t8*(t8*t362+rw*t2*(t499+t2*t3*t321*(t198-t363)))))-t3*t5*t87*(rw*t7*(t3*t10*t190*t534+Itzz*t3*t8*t228*t243*t321)-rw*t2*(t3*t10*t190*t548+t3*t8*t87*t228*t272*t321)+Itzz*Iwyy*t3*t8*t74*t190));
        double et61 = t10*t87*(Itzz*t10*(t609-t1028+t1101+rw*t7*t8*(t369-t716))-Itzz*t3*t5*(t8*t469+rw*t7*(t754+rw*t3*t7*t700)+rw*t2*t775-rw*t2*t3*(t369-t716)))+Itzz*t3*t5*(rw*t2*t1172-rw*t7*t1170);
        double et62 = t1409;
        double et63 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et64 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et65 = t1410;
        double et66 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et67 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et68 = Itzz*(Iwyy*(-t1060+t1141+(t947*(t275-t415))/2.0+rw*t2*t8*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0))+rw*t3*t7*(t836-t1069+t190*t243*(t911+(t611*(t198-t363))/2.0)))+rw*t2*t1350-Iwyy*t8*(t961+t236*(t443+t673)+(t680*(t441+t527))/4.0+rw*t3*t7*(t812+(t518*(t275-t415))/2.0+t190*t272*(t443+t673))))-rw*t8*(Itzz*(t8*t321*(t815+t938+t1200-(t495*(t275-t415))/2.0)+t3*t10*t190*t1339)+Itzz*t3*t5*(rw*t2*t991+t3*t10*t190*t588+Itzz*Iwyy*t3*t8*t228*t243*t321));
        double et69 = Itzz*t10*(Itzz*t10*(t815+t938+t1200-(t495*(t275-t415))/2.0)+t3*t5*t87*t1168)+Itzz*t3*t5*(Iwyy*t1170+rw*t2*t1208);
        double et70 = t1409;
        double et71 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et72 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et73 = rw*t376*(Itzz*(Iwyy*t1362+Iwyy*t8*t1297-rw*t7*t1350)-t10*t87*(Itzz*t10*(-t817+t1198+(t516*(t275-t415))/2.0+t190*t272*(t347-t847))-t3*t5*t87*(t551+rw*t7*t942-rw*t2*t3*t802+(Iwyy*t190*t272*(t276-t491))/2.0))-rw*t8*(Itzz*(t8*t321*(-t817+t1198+(t516*(t275-t415))/2.0+t190*t272*(t347-t847))-t3*t10*t190*(-Iwyy*(t911+(t611*(t198-t363))/2.0)+Iwyy*t8*t757+rw*t7*t1247+rw*t7*t8*t870))-Itzz*t3*t5*(rw*t7*t991+t3*t10*t190*t596+Iwyy*t3*t8*t87*t228*t272*t321))+t3*t5*t87*(Iwyy*t1172+rw*t7*t1208));
        double et74 = 1.0/(et70+et71+et72);
        double et75 = (t997*(Itzz*t1406-t10*t87*t1342+rw*t72*t321*(t426+Itzz*(t347+t529-t772-t847-t950+t1049))+Itzz*t3*t5*(Itzz*t10*t1175-t3*t5*t87*t1173)))/(et56+et57+et58)+t626*t1421*(et59+et60+et61)+t1104*t1422*(t1232+t1396-rw*t3*t8*t10*t190*(t426+Itzz*(t347+t529-t772-t847-t950+t1049)))+(t835*(t10*t1396+rw*t8*(Itzz*t3*t5*t1301+t3*t74*t190*t1287)+Itzz*t3*t5*t1382))/(et62+et63+et64)+((Itzz*t1382+rw*t8*(Itzz*t1301+t30*t59*t74*t228*t370)+t10*t87*t1230)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et65+et66+et67)+rw*t368*t1421*(et68+et69);
        double et76 = et73*et74+rw*t3*t122*t1421*(t8*t321*(t1232+t1396)+t3*t10*t190*(Itzz*t1406-t10*t87*t1342+Itzz*t3*t5*(Itzz*t10*t1175-t3*t5*t87*t1173)));
        double et77 = t1410;
        double et78 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et79 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);

        double et80 = t1410;
        double et81 = t29*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et82 = -t10*t89*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-t3*t5*t29*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et83 = t1409;
        double et84 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et85 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et86 = t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))-rw*t8*(t1384+t8*t321*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))));
        double et87 = rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))));
        double et88 = t1409;
        double et89 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et90 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et91 = t10*t87*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+Itzz*t3*t5*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))));
        double et92 = -Itzz*rw*t3*t5*t8*(t882+Iwyy*(t783-t862)+rw*t2*(t1117+t8*t321*(-t719+t778+t190*t272*(t475-t567))));
        double et93 = t1409;
        double et94 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et95 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et96 = t1409;
        double et97 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et98 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et99 = t1409;
        double et100 = Itzz*(t1326-t1369+t1414+rw*t2*(t1227-t1305+t2*t126*(-t1093+t1148+(t985*(t275-t415))/2.0))+rw*t7*t8*(t1078-t1212+rw*t2*(-t1217+t1225+t238*(-t984+t236*(t499+t2*t3*t321*(t198-t363))+t190*t243*(t636+(t680*(t198-t363))/4.0)))+rw*t3*t7*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+rw*t2*t3*(t1127-t1374+Iwyy*(-t1088+t1096+(t751*(t275-t415))/2.0)+rw*t2*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624)))));
        double et101 = -t10*t87*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))-Itzz*t3*t5*(Itzz*t10*t1346-t3*t5*t87*t1344);
        double et102 = ((t1386+Itzz*rw*t3*t5*t8*t1139)*(-t329-t417+t579+t686+t957+dne*(t472+dom*(t258-t317)+dth*(t301-t352))))/(et77+et78+et79);
        double et103 = (t1104*(Itzz*t10*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))+Itzz*t3*t5*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567))))-Itzz*t8*t10*t59*t238*t1076))/(et80+et81+et82)+(t835*(et86+et87))/(et83+et84+et85);
        double et104 = (t626*(-Itzz*t10*t1381+t3*t5*t87*(t1001-t1222-t1322+rw*t2*t3*(t650+rw*t2*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+Itzz*rw*t3*t5*t8*(-t602+t1068+rw*t7*(t783-t862))))/(et88+et89+et90)+t997*t1422*(t10*t87*(-t896+t916-t1189+t1329+t236*(t347-t847)+rw*t2*(-t1080+t238*(t911+(t611*(t198-t363))/2.0)+t302*(t489-t624)))+Itzz*t3*t5*(-t725+t917+t988+t1165+t1188+t1285)+Itzz*rw*t3*t5*t72*t321*t1076)+(rw*t368*(et91+et92))/(et93+et94+et95);
        double et105 = (rw*t376*(t10*t87*(t1129-t1276-t1375+rw*t7*t8*(t971+Iwyy*(-t809+(t399*(t499+t2*t3*t321*(t198-t363)))/2.0+t190*t243*(t489-t624))))+Itzz*t3*t5*t1394+Itzz*rw*t3*t5*t8*t1337))/(et96+et97+et98);
        double et106 = (rw*t3*t122*(Itzz*t10*(t1384+t8*t321*(t792+t935-t1238+t1318+(t399*(t347-t847))/2.0+rw*t7*t8*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))+Itzz*t3*t5*(t1364+t8*t321*(-t691-t936+rw*t7*(t784-t900+rw*t3*t7*(-t719+t778+t190*t272*(t475-t567)))+(Iwyy*t399*(t276-t491))/4.0+rw*t2*(t806-t899)+rw*t2*t3*(-t606+t717+rw*t2*(-t719+t778+t190*t272*(t475-t567)))))))/(et99+et100+et101);

        d_dNE[index] = et102+et103+et104+et105+et106;
        
        index = index + num_threads;
        // index1 = (grid_size[0] > 1) ? index : 0;
        index2 = (grid_size[1] > 1) ? index : 0;
        index3 = (grid_size[2] > 1) ? index : 0;
        index4 = (grid_size[3] > 1) ? index : 0;
        index5 = (grid_size[4] > 1) ? index : 0;
        index6 = (grid_size[5] > 1) ? index : 0;
        index7 = (grid_size[6] > 1) ? index : 0;
        index8 = (grid_size[7] > 1) ? index : 0;

        uindex1 = (active_actions[0] == 1) ? index : 0;
        uindex2 = (active_actions[1] == 1) ? index : 0;
    }
}

__global__ void step(double* const x1n, // outputs
                     double* const x1, // current states
                     double* const d_x1, // current derivatives
                     const double delta, int32_t const* const grid_size)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3]
                        * grid_size[4] * grid_size[5] * grid_size[6] * grid_size[7];
    
    while (index < num_elements)
    {
        x1n[index] = x1[index] + delta * d_x1[index];
        index = index + num_threads;
    }
}

__global__ void rk4_wbounds(double* const x1n, // outputs
                    double* const x1, // current states
                    double* const k1_x1, // K1
                    double* const k2_x1, // K2
                    double* const k3_x1, // K3
                    double* const k4_x1, // K4
                    double dt, double lower_limit, double upper_limit, int32_t const* const grid_size)
{
    int index = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.x + threadIdx.y * blockDim.x);
    int num_threads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
    int num_elements = grid_size[0] * grid_size[1] * grid_size[2] * grid_size[3]
                        * grid_size[4] * grid_size[5] * grid_size[6] * grid_size[7];
    while (index < num_elements)
    {
        x1n[index] = x1[index] + dt / 6.0 * (k1_x1[index] + 2.0 * k2_x1[index] + 2.0 * k3_x1[index] + k4_x1[index]);        
        x1n[index] = max(lower_limit, min(upper_limit, x1n[index]));
        index = index + num_threads;
    }
}
/*
 * Host code
 */
void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, mxArray const* prhs[])
{
    /* Declare all variables*/
    // Inputs
    // GPU Arrays
    mxGPUArray const* PH;
    mxGPUArray const* TH;
    mxGPUArray const* OM;
    mxGPUArray const* dPH;
    mxGPUArray const* dTH;
    mxGPUArray const* dPS;
    mxGPUArray const* dOM;
    mxGPUArray const* dNE;
    mxGPUArray const* in1;
    mxGPUArray const* in2;
    mxGPUArray const* grid_size;
    mxGPUArray const* active_actions;
    // Pointers for GPU Arrays
    double* p_PH; 
    double* p_TH; 
    double* p_OM;
    double* p_dPH;
    double* p_dTH; 
    double* p_dPS;
    double* p_dOM;
    double* p_dNE;
    double const* p_in1; 
    double const* p_in2; 
    int32_t const* p_grid_size;
    int32_t const* p_active_actions;
    // Pointers for normal inputs
    double* p_mw;
    double* p_Iw;
    double* p_mf;
    double* p_If;
    double* p_mt;
    double* p_It;
    double* p_nw;
    double* p_nt;
    double* p_rw;
    double* p_rf;
    double* p_rt;
    double* p_fcoeff;
    double* p_g;
    double* p_dt;
    double* p_limits;
    double* p_x_dims_free;
    
    // Intermediate variables
    mxGPUArray* k1_PH;
    mxGPUArray* k1_TH;
    mxGPUArray* k1_OM;
    mxGPUArray* k1_dPH;
    mxGPUArray* k1_dTH;
    mxGPUArray* k1_dPS;
    mxGPUArray* k1_dOM;
    mxGPUArray* k1_dNE;
    mxGPUArray* k2_PH;
    mxGPUArray* k2_TH;
    mxGPUArray* k2_OM;
    mxGPUArray* k2_dPH;
    mxGPUArray* k2_dTH;
    mxGPUArray* k2_dPS;
    mxGPUArray* k2_dOM;
    mxGPUArray* k2_dNE;
    mxGPUArray* k3_PH;
    mxGPUArray* k3_TH;
    mxGPUArray* k3_OM;
    mxGPUArray* k3_dPH;
    mxGPUArray* k3_dTH;
    mxGPUArray* k3_dPS;
    mxGPUArray* k3_dOM;
    mxGPUArray* k3_dNE;
    mxGPUArray* k4_PH;
    mxGPUArray* k4_TH;
    mxGPUArray* k4_OM;
    mxGPUArray* k4_dPH;
    mxGPUArray* k4_dTH;
    mxGPUArray* k4_dPS;
    mxGPUArray* k4_dOM;
    mxGPUArray* k4_dNE;

    // Pointers for intermediate variables
    double* p_k1_PH;
    double* p_k1_TH;
    double* p_k1_OM;
    double* p_k1_dPH;
    double* p_k1_dTH;
    double* p_k1_dPS;
    double* p_k1_dOM;
    double* p_k1_dNE;
    double* p_k2_PH;
    double* p_k2_TH;
    double* p_k2_OM;
    double* p_k2_dPH;
    double* p_k2_dTH;
    double* p_k2_dPS;
    double* p_k2_dOM;
    double* p_k2_dNE;
    double* p_k3_PH;
    double* p_k3_TH;
    double* p_k3_OM;
    double* p_k3_dPH;
    double* p_k3_dTH;
    double* p_k3_dPS;
    double* p_k3_dOM;
    double* p_k3_dNE;
    double* p_k4_PH;
    double* p_k4_TH;
    double* p_k4_OM;
    double* p_k4_dPH;
    double* p_k4_dTH;
    double* p_k4_dPS;
    double* p_k4_dOM;
    double* p_k4_dNE;

    // Outputs
    mxGPUArray* PHn;
    mxGPUArray* THn;
    mxGPUArray* OMn;
    mxGPUArray* dPHn;
    mxGPUArray* dTHn;
    mxGPUArray* dPSn;
    mxGPUArray* dOMn;
    mxGPUArray* dNEn;
    // Pointers for outputs
    double* p_PHn;
    double* p_THn;
    double* p_OMn;
    double* p_dPHn;
    double* p_dTHn;
    double* p_dPSn;
    double* p_dOMn;
    double* p_dNEn;

    char const* const errId = "parallel:gpu:mexGPUExample:InvalidInput";
    char const* const errMsg = "Invalid input to MEX file.";

    /* Choose a reasonably sized number of threads for the block. */
    int const threadsPerBlock = 256;
    int const blocksPerGrid = 1024;

    /* Initialize the MathWorks GPU API. */
    mxInitGPU();

    /* Throw an error if the input is not a GPU array. */
    if ((nrhs != 34) || !(mxIsGPUArray(prhs[0])) || !(mxIsGPUArray(prhs[1])) || !(mxIsGPUArray(prhs[2])) || !(mxIsGPUArray(prhs[3])) 
                     || !(mxIsGPUArray(prhs[4])) || !(mxIsGPUArray(prhs[5])) || !(mxIsGPUArray(prhs[6])) || !(mxIsGPUArray(prhs[7]))
                     || !(mxIsGPUArray(prhs[8])) || !(mxIsGPUArray(prhs[9]))) {
        mexErrMsgIdAndTxt(errId, errMsg);
    }

    PH = mxGPUCreateFromMxArray(prhs[0]);
    TH = mxGPUCreateFromMxArray(prhs[1]);
    OM = mxGPUCreateFromMxArray(prhs[2]);
    dPH = mxGPUCreateFromMxArray(prhs[3]);
    dTH = mxGPUCreateFromMxArray(prhs[4]);
    dPS = mxGPUCreateFromMxArray(prhs[5]);
    dOM = mxGPUCreateFromMxArray(prhs[6]);
    dNE = mxGPUCreateFromMxArray(prhs[7]);
    in1 = mxGPUCreateFromMxArray(prhs[8]);
    in2 = mxGPUCreateFromMxArray(prhs[9]);
    grid_size = mxGPUCreateFromMxArray(prhs[30]);
    active_actions = mxGPUCreateFromMxArray(prhs[31]);
    
    p_mw = mxGetDoubles(prhs[10]); 
    double const mw = p_mw[0];
    p_Iw = mxGetDoubles(prhs[11]); 
    double const Iwxx = p_Iw[0];
    p_Iw = mxGetDoubles(prhs[12]); 
    double const Iwyy = p_Iw[0];
    p_Iw = mxGetDoubles(prhs[13]); 
    double const Iwzz = p_Iw[0];

    p_mf = mxGetDoubles(prhs[14]); 
    double const mf = p_mf[0];
    p_If = mxGetDoubles(prhs[15]); 
    double const Ifxx = p_If[0];
    p_If = mxGetDoubles(prhs[16]); 
    double const Ifyy = p_If[0];
    p_If = mxGetDoubles(prhs[17]); 
    double const Ifzz = p_If[0];

    p_mt = mxGetDoubles(prhs[18]); 
    double const mt = p_mt[0];
    p_It = mxGetDoubles(prhs[19]); 
    double const Itxx = p_It[0];
    p_It = mxGetDoubles(prhs[20]); 
    double const Ityy = p_It[0];
    p_It = mxGetDoubles(prhs[21]); 
    double const Itzz = p_It[0];

    p_nw = mxGetDoubles(prhs[22]); 
    double const nw = p_nw[0];
    p_nt = mxGetDoubles(prhs[23]); 
    double const nt = p_nt[0];

    p_rw = mxGetDoubles(prhs[24]); 
    double const rw = p_rw[0];
    p_rf = mxGetDoubles(prhs[25]); 
    double const rf = p_rf[0];
    p_rt = mxGetDoubles(prhs[26]); 
    double const rt = p_rt[0];

    p_fcoeff = mxGetDoubles(prhs[27]);   
    double const fcoeff = p_fcoeff[0];
    p_g = mxGetDoubles(prhs[28]);
    double const g = p_g[0];
    
    p_dt = mxGetDoubles(prhs[29]);
    double const dt = p_dt[0];
    
    p_limits = mxGetDoubles(prhs[32]);
    
    p_x_dims_free = mxGetDoubles(prhs[33]);
    mwSize const* num_x_dims = mxGetDimensions(prhs[33]);
    if ((mxGetNumberOfDimensions(prhs[33]) != 2) || (num_x_dims[1] > 1)) {
        mexErrMsgIdAndTxt(errId, errMsg);
    }
    
    /*
     * Verify that inputs are of appropriate type before extracting the pointer.
     */
    if ((mxGPUGetClassID(PH) != mxDOUBLE_CLASS) || (mxGPUGetClassID(TH) != mxDOUBLE_CLASS) 
        || (mxGPUGetClassID(OM) != mxDOUBLE_CLASS) || (mxGPUGetClassID(dPH) != mxDOUBLE_CLASS) 
        || (mxGPUGetClassID(dTH) != mxDOUBLE_CLASS) || (mxGPUGetClassID(dPS) != mxDOUBLE_CLASS) 
        || (mxGPUGetClassID(dOM) != mxDOUBLE_CLASS) || (mxGPUGetClassID(dNE) != mxDOUBLE_CLASS) 
        || (mxGPUGetClassID(in1) != mxDOUBLE_CLASS) || (mxGPUGetClassID(in2) != mxDOUBLE_CLASS)
        || (mxGPUGetClassID(grid_size) != mxINT32_CLASS) || (mxGPUGetClassID(active_actions) != mxINT32_CLASS)) {
        mexErrMsgIdAndTxt(errId, errMsg);
    }

    /*
     * Now that we have verified the data type, extract a pointer to the input
     * data on the device.
     */
    p_PH = (double*)(mxGPUGetDataReadOnly(PH));
    p_TH = (double*)(mxGPUGetDataReadOnly(TH));
    p_OM = (double*)(mxGPUGetDataReadOnly(OM));
    p_dPH = (double*)(mxGPUGetDataReadOnly(dPH));
    p_dTH = (double*)(mxGPUGetDataReadOnly(dTH));
    p_dPS = (double*)(mxGPUGetDataReadOnly(dPS));
    p_dOM = (double*)(mxGPUGetDataReadOnly(dOM));
    p_dNE = (double*)(mxGPUGetDataReadOnly(dNE));
    p_in1 = (double const*)(mxGPUGetDataReadOnly(in1));
    p_in2 = (double const*)(mxGPUGetDataReadOnly(in2));
    p_grid_size = (int32_t const*)(mxGPUGetDataReadOnly(grid_size));
    p_active_actions = (int32_t const*)(mxGPUGetDataReadOnly(active_actions));
    
    /* Create output arrays*/
    PHn = mxGPUCopyGPUArray(PH);
    p_PHn = (double*)(mxGPUGetData(PHn));

    THn = mxGPUCopyGPUArray(TH);
    p_THn = (double*)(mxGPUGetData(THn));
    
    OMn = mxGPUCopyGPUArray(OM);
    p_OMn = (double*)(mxGPUGetData(OMn));
    
    dPHn = mxGPUCopyGPUArray(dPH);
    p_dPHn = (double*)(mxGPUGetData(dPHn));
    
    dTHn = mxGPUCopyGPUArray(dTH);
    p_dTHn = (double*)(mxGPUGetData(dTHn));
    
    dPSn = mxGPUCopyGPUArray(dPS);
    p_dPSn = (double*)(mxGPUGetData(dPSn));
    
    dOMn = mxGPUCopyGPUArray(dOM);
    p_dOMn = (double*)(mxGPUGetData(dOMn));
    
    dNEn = mxGPUCopyGPUArray(dNE);
    p_dNEn = (double*)(mxGPUGetData(dNEn));
    
    // RK4 - Step 1
    int32_t curr_free_dim;
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            p_k1_PH = p_dPH;
            step << <blocksPerGrid, threadsPerBlock >> > (p_PHn, p_PH, p_k1_PH, 0.5 * dt, p_grid_size);
        }
        else if (curr_free_dim == 2) {
            p_k1_TH = p_dTH;
            step << <blocksPerGrid, threadsPerBlock >> > (p_THn, p_TH, p_k1_TH, 0.5 * dt, p_grid_size);
    
        } else if (curr_free_dim == 3) {
            p_k1_OM = p_dOM;
            step << <blocksPerGrid, threadsPerBlock >> > (p_OMn, p_OM, p_k1_OM, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 4) {
            k1_dPH = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dPH), mxGPUGetDimensions(dPH), mxGPUGetClassID(dPH), mxGPUGetComplexity(dPH), MX_GPU_DO_NOT_INITIALIZE);
            p_k1_dPH = (double*)(mxGPUGetData(k1_dPH));
            dyn3_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k1_dPH, p_PH, p_TH, p_OM, p_dPH, p_dTH, p_dPS, p_dOM, p_dNE, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);
            step << <blocksPerGrid, threadsPerBlock >> > (p_dPHn, p_dPH, p_k1_dPH, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 5) {
            k1_dTH = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dTH), mxGPUGetDimensions(dTH), mxGPUGetClassID(dTH), mxGPUGetComplexity(dTH), MX_GPU_DO_NOT_INITIALIZE);
            p_k1_dTH = (double*)(mxGPUGetData(k1_dTH));
            dyn4_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k1_dTH, p_PH, p_TH, p_OM, p_dPH, p_dTH, p_dPS, p_dOM, p_dNE, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);
            step << <blocksPerGrid, threadsPerBlock >> > (p_dTHn, p_dTH, p_k1_dTH, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 6) {
            k1_dPS = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dPS), mxGPUGetDimensions(dPS), mxGPUGetClassID(dPS), mxGPUGetComplexity(dPS), MX_GPU_DO_NOT_INITIALIZE);
            p_k1_dPS = (double*)(mxGPUGetData(k1_dPS));
            dyn5_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k1_dPS, p_PH, p_TH, p_OM, p_dPH, p_dTH, p_dPS, p_dOM, p_dNE, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);
            step << <blocksPerGrid, threadsPerBlock >> > (p_dPSn, p_dPS, p_k1_dPS, 0.5 * dt, p_grid_size);

        } else if (curr_free_dim == 7) {
            k1_dOM = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dOM), mxGPUGetDimensions(dOM), mxGPUGetClassID(dOM), mxGPUGetComplexity(dOM), MX_GPU_DO_NOT_INITIALIZE);
            p_k1_dOM = (double*)(mxGPUGetData(k1_dOM));
            dyn6_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k1_dOM, p_PH, p_TH, p_OM, p_dPH, p_dTH, p_dPS, p_dOM, p_dNE, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);
            step << <blocksPerGrid, threadsPerBlock >> > (p_dOMn, p_dOM, p_k1_dOM, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 8) {
            k1_dNE = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dNE), mxGPUGetDimensions(dNE), mxGPUGetClassID(dNE), mxGPUGetComplexity(dNE), MX_GPU_DO_NOT_INITIALIZE);
            p_k1_dNE = (double*)(mxGPUGetData(k1_dNE));
            dyn7_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k1_dNE, p_PH, p_TH, p_OM, p_dPH, p_dTH, p_dPS, p_dOM, p_dNE, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);
            step << <blocksPerGrid, threadsPerBlock >> > (p_dNEn, p_dNE, p_k1_dNE, 0.5 * dt, p_grid_size);    
        } 
    }
    
    // RK4 - Step 2
    // Continuous Dynamics
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            k2_PH = mxGPUCopyGPUArray(dPHn);
            p_k2_PH = (double*)(mxGPUGetData(k2_PH));
    
        } else if (curr_free_dim == 2) {
            k2_TH = mxGPUCopyGPUArray(dTHn);
            p_k2_TH = (double*)(mxGPUGetData(k2_TH));
    
        } else if (curr_free_dim == 3) {
            k2_OM = mxGPUCopyGPUArray(dOMn);
            p_k2_OM = (double*)(mxGPUGetData(k2_OM));
            
        } else if (curr_free_dim == 4) {
            k2_dPH = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dPH), mxGPUGetDimensions(dPH), mxGPUGetClassID(dPH), mxGPUGetComplexity(dPH), MX_GPU_DO_NOT_INITIALIZE);
            p_k2_dPH = (double*)(mxGPUGetData(k2_dPH));
            dyn3_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k2_dPH, p_PHn, p_THn, p_OMn, p_dPHn, p_dTHn, p_dPSn, p_dOMn, p_dNEn, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);
            
        } else if (curr_free_dim == 5) {
            k2_dTH = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dTH), mxGPUGetDimensions(dTH), mxGPUGetClassID(dTH), mxGPUGetComplexity(dTH), MX_GPU_DO_NOT_INITIALIZE);
            p_k2_dTH = (double*)(mxGPUGetData(k2_dTH));
            dyn4_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k2_dTH, p_PHn, p_THn, p_OMn, p_dPHn, p_dTHn, p_dPSn, p_dOMn, p_dNEn, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);
            
        } else if (curr_free_dim == 6) {
            k2_dPS = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dPS), mxGPUGetDimensions(dPS), mxGPUGetClassID(dPS), mxGPUGetComplexity(dPS), MX_GPU_DO_NOT_INITIALIZE);
            p_k2_dPS = (double*)(mxGPUGetData(k2_dPS));
            dyn5_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k2_dPS, p_PHn, p_THn, p_OMn, p_dPHn, p_dTHn, p_dPSn, p_dOMn, p_dNEn, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);

        } else if (curr_free_dim == 7) {
            k2_dOM = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dOM), mxGPUGetDimensions(dOM), mxGPUGetClassID(dOM), mxGPUGetComplexity(dOM), MX_GPU_DO_NOT_INITIALIZE);
            p_k2_dOM = (double*)(mxGPUGetData(k2_dOM));
            dyn6_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k2_dOM, p_PHn, p_THn, p_OMn, p_dPHn, p_dTHn, p_dPSn, p_dOMn, p_dNEn, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);
            
        } else if (curr_free_dim == 8) {
            k2_dNE = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dNE), mxGPUGetDimensions(dNE), mxGPUGetClassID(dNE), mxGPUGetComplexity(dNE), MX_GPU_DO_NOT_INITIALIZE);
            p_k2_dNE = (double*)(mxGPUGetData(k2_dNE));
            dyn7_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k2_dNE, p_PHn, p_THn, p_OMn, p_dPHn, p_dTHn, p_dPSn, p_dOMn, p_dNEn, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);    
        } 
    }
    // Discrete step
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_PHn, p_PH, p_k2_PH, 0.5 * dt, p_grid_size);
    
        } else if (curr_free_dim == 2) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_THn, p_TH, p_k2_TH, 0.5 * dt, p_grid_size);
    
        } else if (curr_free_dim == 3) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_OMn, p_OM, p_k2_OM, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 4) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_dPHn, p_dPH, p_k2_dPH, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 5) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_dTHn, p_dTH, p_k2_dTH, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 6) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_dPSn, p_dPS, p_k2_dPS, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 7) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_dOMn, p_dOM, p_k2_dOM, 0.5 * dt, p_grid_size);
            
        } else if (curr_free_dim == 8) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_dNEn, p_dNE, p_k2_dNE, 0.5 * dt, p_grid_size);
        
        } 
    }
    
    // RK4 - Step 3
    // Continuous Dynamics
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            k3_PH = mxGPUCopyGPUArray(dPHn);
            p_k3_PH = (double*)(mxGPUGetData(k3_PH));
    
        } else if (curr_free_dim == 2) {
            k3_TH = mxGPUCopyGPUArray(dTHn);
            p_k3_TH = (double*)(mxGPUGetData(k3_TH));
    
        } else if (curr_free_dim == 3) {
            k3_OM = mxGPUCopyGPUArray(dOMn);
            p_k3_OM = (double*)(mxGPUGetData(k3_OM));
            
        } else if (curr_free_dim == 4) {
            k3_dPH = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dPH), mxGPUGetDimensions(dPH), mxGPUGetClassID(dPH), mxGPUGetComplexity(dPH), MX_GPU_DO_NOT_INITIALIZE);
            p_k3_dPH = (double*)(mxGPUGetData(k3_dPH));
            dyn3_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k3_dPH, p_PHn, p_THn, p_OMn, p_dPHn, p_dTHn, p_dPSn, p_dOMn, p_dNEn, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);
            
        } else if (curr_free_dim == 5) {
            k3_dTH = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dTH), mxGPUGetDimensions(dTH), mxGPUGetClassID(dTH), mxGPUGetComplexity(dTH), MX_GPU_DO_NOT_INITIALIZE);
            p_k3_dTH = (double*)(mxGPUGetData(k3_dTH));
            dyn4_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k3_dTH, p_PHn, p_THn, p_OMn, p_dPHn, p_dTHn, p_dPSn, p_dOMn, p_dNEn, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);
            
        } else if (curr_free_dim == 6) {
            k3_dPS = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dPS), mxGPUGetDimensions(dPS), mxGPUGetClassID(dPS), mxGPUGetComplexity(dPS), MX_GPU_DO_NOT_INITIALIZE);
            p_k3_dPS = (double*)(mxGPUGetData(k3_dPS));
            dyn5_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k3_dPS, p_PHn, p_THn, p_OMn, p_dPHn, p_dTHn, p_dPSn, p_dOMn, p_dNEn, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);

        } else if (curr_free_dim == 7) {
            k3_dOM = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dOM), mxGPUGetDimensions(dOM), mxGPUGetClassID(dOM), mxGPUGetComplexity(dOM), MX_GPU_DO_NOT_INITIALIZE);
            p_k3_dOM = (double*)(mxGPUGetData(k3_dOM));
            dyn6_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k3_dOM, p_PHn, p_THn, p_OMn, p_dPHn, p_dTHn, p_dPSn, p_dOMn, p_dNEn, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);
            
        } else if (curr_free_dim == 8) {
            k3_dNE = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dNE), mxGPUGetDimensions(dNE), mxGPUGetClassID(dNE), mxGPUGetComplexity(dNE), MX_GPU_DO_NOT_INITIALIZE);
            p_k3_dNE = (double*)(mxGPUGetData(k3_dNE));
            dyn7_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k3_dNE, p_PHn, p_THn, p_OMn, p_dPHn, p_dTHn, p_dPSn, p_dOMn, p_dNEn, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);    
        } 
    }
    // Discrete step
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_PHn, p_PH, p_k3_PH, dt, p_grid_size);

        } else if (curr_free_dim == 2) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_THn, p_TH, p_k3_TH, dt, p_grid_size);
    
        } else if (curr_free_dim == 3) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_OMn, p_OM, p_k3_OM, dt, p_grid_size);
            
        } else if (curr_free_dim == 4) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_dPHn, p_dPH, p_k3_dPH, dt, p_grid_size);
            
        } else if (curr_free_dim == 5) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_dTHn, p_dTH, p_k3_dTH, dt, p_grid_size);
            
        } else if (curr_free_dim == 6) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_dPSn, p_dPS, p_k3_dPS, dt, p_grid_size);
            
        } else if (curr_free_dim == 7) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_dOMn, p_dOM, p_k3_dOM, dt, p_grid_size);
            
        } else if (curr_free_dim == 8) {
            step << <blocksPerGrid, threadsPerBlock >> > (p_dNEn, p_dNE, p_k3_dNE, dt, p_grid_size);
        
        } 
    }
    
    // RK4 - Step 4
    // Continuous Dynamics
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            k4_PH = mxGPUCopyGPUArray(dPHn);
            p_k4_PH = (double*)(mxGPUGetData(k4_PH));
    
        } else if (curr_free_dim == 2) {
            k4_TH = mxGPUCopyGPUArray(dTHn);
            p_k4_TH = (double*)(mxGPUGetData(k4_TH));
    
        } else if (curr_free_dim == 3) {
            k4_OM = mxGPUCopyGPUArray(dOMn);
            p_k4_OM = (double*)(mxGPUGetData(k4_OM));
            
        } else if (curr_free_dim == 4) {
            k4_dPH = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dPH), mxGPUGetDimensions(dPH), mxGPUGetClassID(dPH), mxGPUGetComplexity(dPH), MX_GPU_DO_NOT_INITIALIZE);
            p_k4_dPH = (double*)(mxGPUGetData(k4_dPH));
            dyn3_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k4_dPH, p_PHn, p_THn, p_OMn, p_dPHn, p_dTHn, p_dPSn, p_dOMn, p_dNEn, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);

        } else if (curr_free_dim == 5) {
            k4_dTH = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dTH), mxGPUGetDimensions(dTH), mxGPUGetClassID(dTH), mxGPUGetComplexity(dTH), MX_GPU_DO_NOT_INITIALIZE);
            p_k4_dTH = (double*)(mxGPUGetData(k4_dTH));
            dyn4_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k4_dTH, p_PHn, p_THn, p_OMn, p_dPHn, p_dTHn, p_dPSn, p_dOMn, p_dNEn, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);
            
        } else if (curr_free_dim == 6) {
            k4_dPS = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dPS), mxGPUGetDimensions(dPS), mxGPUGetClassID(dPS), mxGPUGetComplexity(dPS), MX_GPU_DO_NOT_INITIALIZE);
            p_k4_dPS = (double*)(mxGPUGetData(k4_dPS));
            dyn5_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k4_dPS, p_PHn, p_THn, p_OMn, p_dPHn, p_dTHn, p_dPSn, p_dOMn, p_dNEn, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);

        } else if (curr_free_dim == 7) {
            k4_dOM = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dOM), mxGPUGetDimensions(dOM), mxGPUGetClassID(dOM), mxGPUGetComplexity(dOM), MX_GPU_DO_NOT_INITIALIZE);
            p_k4_dOM = (double*)(mxGPUGetData(k4_dOM));
            dyn6_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k4_dOM, p_PHn, p_THn, p_OMn, p_dPHn, p_dTHn, p_dPSn, p_dOMn, p_dNEn, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);
            
        } else if (curr_free_dim == 8) {
            k4_dNE = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(dNE), mxGPUGetDimensions(dNE), mxGPUGetClassID(dNE), mxGPUGetComplexity(dNE), MX_GPU_DO_NOT_INITIALIZE);
            p_k4_dNE = (double*)(mxGPUGetData(k4_dNE));
            dyn7_mex_continuous << <blocksPerGrid, threadsPerBlock >> > (p_k4_dNE, p_PHn, p_THn, p_OMn, p_dPHn, p_dTHn, p_dPSn, p_dOMn, p_dNEn, p_in1, p_in2,
                                                                         mw, Iwxx, Iwyy, Iwzz, mf, Ifxx, Ifyy, Ifzz, mt, Itxx, Ityy, Itzz, 
                                                                         nw, nt, rw, rf, rt, fcoeff, g, p_grid_size, p_active_actions);    
        } 
    }

    // RK4 Ouput
    for (int32_t ii=0; ii < num_x_dims[0]; ii++) {
        curr_free_dim = int32_t(p_x_dims_free[ii]);
        if (curr_free_dim == 1) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_PHn, p_PH, p_k1_PH, p_k2_PH, p_k3_PH, p_k4_PH, dt, p_limits[0], p_limits[8], p_grid_size);
            mxGPUDestroyGPUArray(k1_PH);
            mxGPUDestroyGPUArray(k2_PH);
            mxGPUDestroyGPUArray(k3_PH);
            mxGPUDestroyGPUArray(k4_PH);
        } else if (curr_free_dim == 2) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_THn, p_TH, p_k1_TH, p_k2_TH, p_k3_TH, p_k4_TH, dt, p_limits[1], p_limits[9], p_grid_size);
            mxGPUDestroyGPUArray(k1_TH);
            mxGPUDestroyGPUArray(k2_TH);
            mxGPUDestroyGPUArray(k3_TH);
            mxGPUDestroyGPUArray(k4_TH);
        } else if (curr_free_dim == 3) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_OMn, p_OM, p_k1_OM, p_k2_OM, p_k3_OM, p_k4_OM, dt, p_limits[2], p_limits[10], p_grid_size);
            mxGPUDestroyGPUArray(k1_OM);
            mxGPUDestroyGPUArray(k2_OM);
            mxGPUDestroyGPUArray(k3_OM);
            mxGPUDestroyGPUArray(k4_OM);
        } else if (curr_free_dim == 4) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_dPHn, p_dPH, p_k1_dPH, p_k2_dPH, p_k3_dPH, p_k4_dPH, dt, p_limits[3], p_limits[11], p_grid_size);
            mxGPUDestroyGPUArray(k1_dPH);
            mxGPUDestroyGPUArray(k2_dPH);
            mxGPUDestroyGPUArray(k3_dPH);
            mxGPUDestroyGPUArray(k4_dPH);
        } else if (curr_free_dim == 5) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_dTHn, p_dTH, p_k1_dTH, p_k2_dTH, p_k3_dTH, p_k4_dTH, dt, p_limits[4], p_limits[12], p_grid_size);
            mxGPUDestroyGPUArray(k1_dTH);
            mxGPUDestroyGPUArray(k2_dTH);
            mxGPUDestroyGPUArray(k3_dTH);
            mxGPUDestroyGPUArray(k4_dTH);
        } else if (curr_free_dim == 6) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_dPSn, p_dPS, p_k1_dPS, p_k2_dPS, p_k3_dPS, p_k4_dPS, dt, p_limits[5], p_limits[13], p_grid_size);
            mxGPUDestroyGPUArray(k1_dPS);
            mxGPUDestroyGPUArray(k2_dPS);
            mxGPUDestroyGPUArray(k3_dPS);
            mxGPUDestroyGPUArray(k4_dPS);
        } else if (curr_free_dim == 7) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_dOMn, p_dOM, p_k1_dOM, p_k2_dOM, p_k3_dOM, p_k4_dOM, dt, p_limits[6], p_limits[14], p_grid_size);
            mxGPUDestroyGPUArray(k1_dOM);
            mxGPUDestroyGPUArray(k2_dOM);
            mxGPUDestroyGPUArray(k3_dOM);
            mxGPUDestroyGPUArray(k4_dOM);
        } else if (curr_free_dim == 8) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_dNEn, p_dNE, p_k1_dNE, p_k2_dNE, p_k3_dNE, p_k4_dNE, dt, p_limits[7], p_limits[15], p_grid_size);
            mxGPUDestroyGPUArray(k1_dNE);
            mxGPUDestroyGPUArray(k2_dNE);
            mxGPUDestroyGPUArray(k3_dNE);
            mxGPUDestroyGPUArray(k4_dNE);
        } 
    }
 
    /* Wrap the result up as a MATLAB gpuArray for return. */
    plhs[0] = mxGPUCreateMxArrayOnGPU(PHn);
    plhs[1] = mxGPUCreateMxArrayOnGPU(THn);
    plhs[2] = mxGPUCreateMxArrayOnGPU(OMn);
    plhs[3] = mxGPUCreateMxArrayOnGPU(dPHn);
    plhs[4] = mxGPUCreateMxArrayOnGPU(dTHn);
    plhs[5] = mxGPUCreateMxArrayOnGPU(dPSn);
    plhs[6] = mxGPUCreateMxArrayOnGPU(dOMn);
    plhs[7] = mxGPUCreateMxArrayOnGPU(dNEn);

    /*
     * The mxGPUArray pointers are host-side structures that refer to device
     * data. These must be destroyed before leaving the MEX function.
     */
    mxGPUDestroyGPUArray(PH);
    mxGPUDestroyGPUArray(TH);
    mxGPUDestroyGPUArray(OM);
    mxGPUDestroyGPUArray(dPH);
    mxGPUDestroyGPUArray(dTH);
    mxGPUDestroyGPUArray(dPS);
    mxGPUDestroyGPUArray(dOM);
    mxGPUDestroyGPUArray(dNE);
    mxGPUDestroyGPUArray(in1);
    mxGPUDestroyGPUArray(in2);
    mxGPUDestroyGPUArray(grid_size);
    mxGPUDestroyGPUArray(active_actions);
    mxGPUDestroyGPUArray(PHn);
    mxGPUDestroyGPUArray(THn);
    mxGPUDestroyGPUArray(OMn);
    mxGPUDestroyGPUArray(dPHn);
    mxGPUDestroyGPUArray(dTHn);
    mxGPUDestroyGPUArray(dPSn);
    mxGPUDestroyGPUArray(dOMn);
    mxGPUDestroyGPUArray(dNEn);
}
