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

    int index1 = (grid_size[0] > 1) ? index : 0;
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

        double t2 = cos(th);
        double t3 = cos(ps);
        double t4 = cos(om);
        double t5 = cos(ne);
        double t6 = sin(th);
        double t7 = sin(ps);
        double t8 = sin(om);
        double t9 = sin(ne);
        double t10 = mf+mt;
        double t11 = mf*rf;
        double t12 = mt*rt;
        double t13 = rf+rt;
        double t14 = Ifxx*2.0;
        double t15 = Ifxx*4.0;
        double t16 = Ifyy*4.0;
        double t17 = Ifzz*2.0;
        double t18 = Ifzz*4.0;
        double t19 = Itxx*2.0;
        double t20 = Itxx*3.0;
        double t21 = Itxx*4.0;
        double t22 = Ityy*2.0;
        double t23 = Ityy*3.0;
        double t24 = Ityy*4.0;
        double t25 = Itzz*2.0;
        double t26 = Itzz*4.0;
        double t27 = Itzz*Itzz;
        double t28 = Iwxx*2.0;
        double t29 = Iwxx*4.0;
        double t30 = Iwyy*4.0;
        double t31 = Iwzz*2.0;
        double t32 = Iwzz*4.0;
        double t33 = mf*2.0;
        double t34 = mt*2.0;
        double t35 = mt*4.0;
        double t36 = mw*2.0;
        double t37 = mw*4.0;
        double t38 = rf*rf;

        double t40 = rw*rw;
        double t41 = th*2.0;
        double t42 = ps*2.0;
        double t43 = om*2.0;
        double t44 = ne*2.0;
        double t45 = dps*2.0;
        double t46 = dom*2.0;
        double t47 = dph*dph;
        double t48 = dth*dth;
        double t49 = dom*dom;
        double t73 = Ifxx*8.0;
        double t74 = Ifyy*8.0;
        double t75 = -Ifzz;
        double t78 = Ifzz*8.0;
        double t79 = -Ityy;
        double t83 = -Itzz;
        double t86 = Itzz*8.0;
        double t89 = Iwxx*8.0;
        double t90 = Iwyy*8.0;
        double t91 = -Iwzz;
        double t94 = Iwzz*8.0;
        double t97 = -Tneta;
        double t50 = cos(t41);
        double t51 = cos(t42);
        double t52 = cos(t43);
        double t53 = cos(t44);
        double t54 = t2*t2;
        double t55 = t3*t3;
        double t56 = t4*t4;
        double t57 = t4*t4*t4;
        double t58 = t5*t5;
        double t59 = t11*2.0;
        double t60 = t11*4.0;
        double t61 = t12*2.0;
        double t62 = sin(t41);
        double t63 = sin(t42);
        double t64 = sin(t43);
        double t65 = sin(t44);
        double t66 = t6*t6;
        double t67 = t7*t7;
        double t68 = t8*t8;
        double t69 = t9*t9;
        double t70 = mw+t10;
        double t71 = -t14;
        double t72 = -t15;
        double t76 = -t17;
        double t77 = -t18;
        double t80 = -t22;
        double t81 = -t23;
        double t82 = -t24;
        double t84 = -t25;
        double t85 = -t26;
        double t87 = -t28;
        double t88 = -t29;
        double t92 = -t31;
        double t93 = -t32;
        double t95 = rw*t11;
        double t96 = t6*dph;
        double t98 = rf*t10;
        double t99 = mt*t13;
        double t100 = t5*t6;
        double t101 = t40*2.0;
        double t102 = t6*t9;
        double t103 = t11*8.0;
        double t104 = t12*8.0;
        double t105 = rf*t11;
        double t106 = mt*t38;
        double t107 = mf*t40;
        double t108 = rt*t12;
        double t109 = mt*t40;
        double t110 = mw*t40;
        double t111 = t13*t13;
        double t112 = -t78;
        double t113 = -t86;
        double t114 = -t94;
        double t119 = t13*t34;
        double t120 = t13*t35;
        double t121 = -t48;
        double t123 = Itxx+t79;
        double t124 = Iwxx+t91;
        double t142 = t2*t5*t8;
        double t145 = t2*t8*t9;
        double t152 = t34+t36;
        double t153 = t35+t37;
        double t154 = Ifxx*t2*t4*t8;

        double t157 = Iwxx*t2*t3*t7;

        double t180 = t2*t4*t8*t75;
        double t181 = t2*t4*t8*t83;
        double t182 = t2*t3*t7*t91;
        double t188 = t10*2.0+t36;
        double t115 = rf*t61;
        double t116 = t96*2.0;
        double t117 = t96*4.0;
        double t118 = t98*2.0;
        double t122 = t46*t96;
        double t125 = rw*t70;
        double t126 = rw*t99;
        double t127 = t96+dps;
        double t128 = t96+dom;
        double t129 = t50*3.0;
        double t130 = rf*t59;
        double t131 = t105*3.0;
        double t132 = rf*t60;
        double t133 = -t96;
        double t134 = Itxx*t58;
        double t135 = Ityy*t58;
        double t136 = Iwxx*t55;
        double t137 = Ifzz*t68;
        double t138 = Itxx*t69;
        double t139 = Ityy*t69;
        double t140 = Itzz*t68;
        double t141 = Iwzz*t67;
        double t143 = t98*8.0;
        double t144 = t99*8.0;
        double t146 = t13*t99;
        double t149 = rf*t103;
        double t150 = t16*t66;
        double t151 = t30*t66;
        double t159 = t34*t111;
        double t161 = t35*t111;
        double t162 = t47*t50;
        double t163 = t47*t54;
        double t165 = t19+t80;
        double t166 = t20+t81;
        double t167 = t21+t82;
        double t168 = t123*t123;
        double t169 = t28+t92;
        double t170 = t29+t93;
        double t171 = t40*t70;
        double t172 = t45+t96;
        double t173 = t46+t96;
        double t174 = -t142;
        double t175 = t38+t101;
        double t177 = t12+t98;
        double t178 = t11+t99;
        double t187 = t89+t114;
        double t189 = t52*t54*2.0;
        double t191 = t18*t54*t56;
        double t192 = t26*t54*t56;
        double t193 = t32*t54*t55;
        double t194 = t40*t152;
        double t195 = t40*t153;
        double t196 = t15*t54*t68;
        double t197 = t29*t54*t67;
        double t200 = t59+t119;
        double t201 = t60+t120;
        double t202 = t53*t123;
        double t203 = t51*t124;
        double t204 = t63*t124;
        double t208 = t8*t53*t62*4.0;
        double t219 = t100+t145;
        double t224 = t6*t65*t123;
        double t230 = t5*t9*t49*t123;
        double t231 = t3*t7*t48*t124;
        double t234 = t40*t188;

        double t254 = t2*t4*t65*t123*2.0;
        double t255 = t62*t65*t123*dph;
        double t258 = t4*t62*t65*t123;
        double t259 = t2*t64*t65*t123;

        double t291 = t5*t9*t56*t121*t123;
        double t295 = t2*t5*t9*t56*t83*t123;
        double t147 = t117+dps;
        double t148 = -t129;
        double t160 = t146*3.0;

        double t176 = t13*t144;
        double t179 = t171*4.0;
        double t183 = t96*t127;
        double t184 = t46+t116+dps;
        double t185 = -t162;
        double t186 = Iwyy+t171;
        double t190 = t178*t178;
        double t198 = mf*t175;
        double t199 = t61+t118;
        double t205 = t194*4.0;
        double t206 = t33*t175;
        double t209 = t4*t178;
        double t210 = t202*2.0;
        double t211 = t202*4.0;
        double t212 = t203*4.0;
        double t213 = t53*t166;
        double t214 = rw*t4*t177;
        double t215 = t63*t169;
        double t217 = t105+t146;
        double t218 = Itzz+t202;
        double t220 = t104+t143;
        double t221 = t103+t144;

        double t225 = -t202;
        double t229 = t51*t169*2.0;
        double t233 = t224*4.0;

        double t236 = g*t2*t8*t178;
        double t240 = rw*t4*t200;
        double t241 = rw*t8*t200;
        double t242 = rw*t8*t201;
        double t243 = t130+t159;
        double t244 = rw*t201*dth*dom;
        double t246 = t102+t174;
        double t248 = t4*t224;
        double t250 = t219*t219;
        double t251 = rw*t8*t49*t177;
        double t253 = t54*t204*dps;
        double t257 = t2*t4*t65*t165;
        double t260 = t65*t96*t165*dth;
        double t262 = t9*t123*t174;

        double t264 = t63*t187*dth*dps;
        double t265 = t65*t123*t128;
        double t269 = t56*t65*t165*dth;
        double t270 = t259*2.0;
        double t275 = t3*t7*t124*t163;
        double t276 = t259*dph;
        double t280 = Ityy*t4*t9*t219;
        double t282 = rw*t8*t133*t177*dom;
        double t297 = t8*t50*t65*t167*dph;
        double t306 = t56*t58*t69*t168;
        double t316 = t2*t51*t170*t172*dph;
        double t327 = Itzz*rw*t2*t5*t9*t57*t123*t178;
        double t335 = Ifxx+t105+t106+t108+t115+t134+t139;
        double t336 = Ifyy+t105+t106+t108+t115+t135+t138;
        double t207 = t198*4.0;
        double t216 = rw*t209;
        double t223 = t30+t179;
        double t226 = -t210;

        double t228 = t213*2.0;
        double t232 = t48+t183;
        double t238 = -t213;
        double t245 = Itzz+t225;

        double t252 = t4*t217;
        double t256 = t25+t210;
        double t261 = -t236;

        double t268 = t246*t246;
        double t271 = -t251;
        double t272 = rw*t6*t8*t199*2.0;
        double t273 = fcoeff+t253;
        double t277 = t24*t250;
        double t278 = rw*t2*t220*dth;

        double t284 = t40*t56*t190;
        double t285 = t66*t242;
        double t286 = t2*t8*t218*4.0;
        double t287 = t125+t209;
        double t288 = t148+t189+1.0;
        double t289 = t6*t45*t241;
        double t290 = -t269;
        double t292 = -t275;
        double t296 = t49+t122+t185;
        double t305 = rw*t4*t66*t221;
        double t307 = Itxx*t4*t5*t246;
        double t308 = -t297;

        double t314 = t52*t54*t243;
        double t317 = rw*t2*t4*t147*t201;

        double t328 = rw*t2*t4*t184*t221*dph;
        double t329 = t146+t194+t198;
        double t331 = t186*t295;
        double t337 = t186+t203+t214;
        double t349 = t47*t123*t219*t246;
        double t353 = t56*t335;
        double t354 = Itxx+Ityy+t14+t76+t84+t202+t243;
        double t360 = t75+t83+t335;
        double t365 = t15+t19+t22+t77+t85+t132+t161+t210;
        double t369 = t186*t336;
        double t376 = t21+t24+t73+t112+t113+t149+t176+t211;
        double t239 = -t228;
        double t247 = t223*dps;
        double t267 = t25+t226;
        double t283 = t4*t245*dth;
        double t293 = t21*t268;
        double t294 = t2*t8*t245;
        double t298 = t4*t6*t256;
        double t301 = t273*dph*4.0;

        double t303 = rw*t287;
        double t313 = -t307;
        double t315 = t65*t288;
        double t318 = g*t6*t287*8.0;
        double t320 = t4*t128*t256;
        double t321 = -t314;

        double t323 = rw*t200*t232;
        double t324 = t95+t126+t252;
        double t326 = -dps*(t212-t223);
        double t332 = t65*t123*t296;

        double t339 = t254+t272;
        double t340 = t50*t329;
        double t341 = t255+t278;
        double t343 = t16+t19+t22+t132+t161+t226;

        double t355 = t2*t337*dph*dth;
        double t364 = t64*t354;
        double t371 = t64*t360;
        double t372 = t216+t336;
        double t381 = (t2*t52*t365*dph*dth)/4.0;
        double t382 = t64*t376*dth*dom;

        double t385 = t2*t52*t173*t365*dph;
        double t387 = t2*t52*t365*(t96-dom);
        double t403 = Itxx+Ityy+t16+t30+t71+t76+t84+t87+t92+t159+t195+t206+t238;

        double t415 = t107+t109+t110+t136+t137+t140+t141+t240+t353;
        double t299 = t294*dph;
        double t300 = -t283;
        double t304 = t4*t267*dom;
        double t309 = t4*t6*t267*2.0;
        double t310 = Iwyy+t303;
        double t333 = t2*t8*t324;
        double t334 = t224+t294;
        double t344 = -t340;

        double t350 = t343*dom;
        double t351 = t259+t298;
        double t352 = t49*t339;
        double t361 = (t4*t341*dph)/4.0;
        double t367 = (t2*t343*dph*dth)/4.0;
        double t368 = t260+t323;
        double t370 = -t123*dph*(t208-t315);
        double t373 = t244+t332;
        double t374 = t54*t364;

        double t377 = t6*t372;

        double t390 = (t364*(t48-t163))/4.0;
        double t392 = t276+t290+t320;

        double t401 = t215+t242+t364;
        double t402 = t204+t241+t371;
        double t408 = t96*t403;
        double t419 = t19+t22+t72+t74+t77+t85+t88+t90+t93+t161+t205+t207+t229+t239;
        double t423 = Itzz*t2*t4*t415;
        double t431 = t231+t271+t282+t292+t355+Tneta;
        double t433 = t284*t415;

        double t440 = t336*t415;

        double t319 = t310*t310;
        double t342 = t334*dph*dom;
        double t346 = t4*t9*t100*t123*t310;
        double t347 = t6*t216*t310;
        double t357 = t351*dph;
        double t358 = -t352;

        double t378 = -t374;
        double t379 = (t8*t368)/2.0;
        double t380 = t8*t373*4.0;
        double t384 = t265+t299+t300;
        double t386 = t6*t310*t336;
        double t395 = t392*dne*4.0;
        double t397 = t262+t377;
        double t404 = t6*t401;
        double t405 = t2*t402;
        double t424 = t96*t419;
        double t437 = t154+t157+t180+t181+t182+t280+t313+t333;

        double t448 = t247+t350+t408;
        double t449 = dne*(t370+dom*(t233-t286)+dth*(t270-t309));

        double t481 = t131+t150+t151+t160+t191+t192+t193+t196+t197+t234+t277+t293+t305+t321+t344;

        double t330 = t66*t319;
        double t356 = -t347;

        double t363 = t8*t83*t347;
        double t388 = t384*dne;

        double t394 = t8*t83*t386;
        double t398 = -t395;
        double t399 = dth*(t304-t357)*(-1.0/2.0);
        double t400 = Itzz*t8*t397;
        double t407 = t186*t397;
        double t409 = t4*t5*t9*t123*t397;
        double t410 = t216*t397;
        double t412 = t258+t285+t378;
        double t416 = t257+t404;
        double t417 = t6*t310*t397;

        double t426 = Itzz*t2*t4*(t248-t405)*(-1.0/2.0);
        double t429 = t4*t5*t9*t123*(t248-t405)*(-1.0/2.0);
        double t430 = t216*(t248-t405)*(-1.0/2.0);
        double t439 = Itzz*t8*t437;

        double t442 = t347*t415;
        double t444 = t336*(t248-t405)*(-1.0/2.0);
        double t446 = t4*t5*t9*t123*t437;
        double t447 = t216*t437;
        double t450 = t2*t448*dph*2.0;
        double t457 = t336*t437;
        double t458 = t397*t415;

        double t463 = t326+t350+t424;
        double t475 = t397*t437;

        double t482 = (t437*(t248-t405))/2.0;
        double t485 = (Itzz*t8*t481)/4.0;
        double t487 = (t186*t481)/4.0;
        double t489 = (t4*t5*t9*t123*t481)/4.0;
        double t490 = (t216*t481)/4.0;

        double t512 = (t415*t481)/4.0;
        double t411 = -t410;
        double t413 = t412*dph*2.0;
        double t414 = t295+t400;
        double t421 = t48*t416;

        double t438 = -Itzz*t8*(t347-t407);
        double t445 = t356*t415;
        double t454 = t346+t430;
        double t464 = t2*t463;
        double t466 = t346+t447;
        double t469 = t230+t291+t342+t349+t399+Tomega;

        double t478 = t409+t444;
        double t479 = ((t248-t405)*(t347-t407))/2.0;
        double t484 = t423+t439;
        double t488 = -t487;
        double t491 = -t490;
        double t492 = t409+t457;
        double t498 = t429+t458;
        double t508 = t446+t458;
        double t513 = t426+t485;
        double t516 = t97+t261+t361+t367+t379+t381+t388+t390;

        double t525 = t4*t5*t9*t123*(t489-(t397*(t248-t405))/2.0);
        double t527 = t475+t489;
        double t537 = t264+t316+t318+t328+t380+t382+t385+t398+t450;
        double t539 = t482+t512;

        double t422 = t8*t83*t414;
        double t425 = t186*t414;
        double t427 = t216*t414;
        double t428 = t289+t413;
        double t451 = t386+t411;
        double t456 = t331+t438;
        double t459 = t216*t454;
        double t470 = t4*t5*t9*t123*t466;
        double t471 = t216*t466;
        double t476 = t6*t310*t466;
        double t483 = t186*t478;
        double t486 = Itzz*t2*t4*t484;
        double t495 = t186*t492;
        double t496 = t6*t310*t492;
        double t497 = t330+t488;
        double t499 = t308+t317+t387+t464;
        double t504 = t186*t498;
        double t506 = t216*t498;
        double t509 = t6*t310*t498;
        double t511 = t186*t508;
        double t514 = t417+t491;
        double t515 = Itzz*t8*t513;
        double t520 = t397*t508;
        double t529 = t186*t527;
        double t531 = t216*t527;

        double t541 = t83*t539;
        double t542 = t216*t539;
        double t544 = t336*t539;
        double t432 = t428*dom;
        double t452 = t363+t425;
        double t453 = Itzz*t8*t451;

        double t465 = t394+t427;

        double t493 = t415*t451;

        double t500 = t4*t5*t9*t123*t497;
        double t502 = t499*dth;

        double t510 = t336*t497;
        double t517 = t4*t5*t9*t123*t514;
        double t518 = t216*t514;

        double t533 = -t531;
        double t534 = t445+t511;
        double t543 = -t542;
        double t545 = -t544;

        double t559 = t486+t515+t541;
        double t434 = -t432;
        double t467 = t327+t453;
        double t503 = -t500;

        double t549 = t496+t533;

        double t556 = t509+t543;
        double t560 = t216*t559;

        double t567 = Itzz*t8*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))*8.0;
        double t568 = t520+t525+t545;
        double t538 = t479+t503;
        double t548 = t301+t358+t421+t434+t449+t502;

        double t557 = t216*t556;

        double t569 = t186*t568;

        double t576 = -1.0/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double t577 = -1.0/(t26*(t557-t569+t6*t310*(t470-t493))+t8*t85*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+t2*t4*t26*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double t579 = -1.0/(-t567+t86*(t557-t569+t6*t310*(t470-t493))+t2*t4*t86*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));

        double et1 = (t548*(Itzz*(t433+t186*(t306-t440))-t27*t68*(t284-t369)))/(t26*(t557-t569+t6*t310*(t470-t493))+t8*t85*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+t2*t4*t26*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t431*t576*(Itzz*(t470-t493)-Itzz*t8*t465);
        double et2 = (t469*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t537*t579*(Itzz*(t471-t495)+t2*t4*t8*t27*(t284-t369));
        double et3 = (t516*(Itzz*t534+t8*t83*t452))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et4 = (t537*(Itzz*(-t510+t518+t397*(t347-t407))+t27*t54*t56*(t284-t369)))/(-t567+t86*(t557-t569+t6*t310*(t470-t493))+t2*t4*t86*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+(t516*(Itzz*t538+Itzz*t2*t4*t452))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et5 = (t431*(Itzz*(t517-(t451*(t248-t405))/2.0)+t2*t4*t83*t465))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t469*t576*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483));
        double et6 = (t548*(Itzz*(t459-t483)+t2*t4*t8*t27*(t284-t369)))/(t26*(t557-t569+t6*t310*(t470-t493))+t8*t85*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+t2*t4*t26*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et7 = (t431*(Itzz*t568+t8*t83*(t397*t414-t336*t513)+Itzz*t2*t4*(t336*t484+t4*t5*t9*t123*t414)))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t516*t576*(t560+t6*t310*(t422+Itzz*t508))+t537*t579*(Itzz*t549+Itzz*t2*t4*t467)+t548*t577*(Itzz*(t506+t6*t310*(t306-t440))-t8*t83*t467);
        double et8 = (t469*(Itzz*t8*t549-Itzz*t2*t4*(t506+t6*t310*(t306-t440))))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et9 = (t516*(t186*t559+t330*(Itzz*t107+Itzz*t109+Itzz*t110+Itzz*t136+Itzz*t137+Itzz*t141+Ifxx*Itzz*t56+Itzz*t56*t105+Itzz*t56*t106+Itzz*t56*t108+Itzz*t56*t134+Itzz*t56*t139+t4*t25*t95+rf*t12*t25*t56+rw*t4*t12*t25+mt*rf*rw*t4*t25)))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et10 = (t537*(Itzz*(t476-t529)-Itzz*t2*t4*t456))/(-t567+t86*(t557-t569+t6*t310*(t470-t493))+t2*t4*t86*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t431*t576*(t560+t6*t310*(t422+Itzz*t498))+t469*t576*(Itzz*t8*(t476-t529)+Itzz*t2*t4*(t442-t504))+t548*t577*(Itzz*(t442-t504)-t8*t83*t456);
        double et11 = t516*t576*(Itzz*t8*t538+Itzz*t2*t4*t534)+(t469*(t557-t569+t6*t310*(t470-t493)))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t537*t579*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t471-t495))+t548*t577*(Itzz*t8*(t459-t483)+Itzz*t2*t4*(t433+t186*(t306-t440)));
        double et12 = t431*t576*(Itzz*t8*(t517-(t451*(t248-t405))/2.0)-Itzz*t2*t4*(t470-t493));

        d_dPH[index] = et1+et2+et3;
        
        index = index + num_threads;
        index1 = (grid_size[0] > 1) ? index : 0;
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

    int index1 = (grid_size[0] > 1) ? index : 0;
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
        th = TH[index1];
        om = OM[index2];
        dph = dPH[index3];
        dth = dTH[index4];
        dps = dPS[index5];
        dom = dOM[index6];
        dne = dNE[index7];

        Tomega = nw*in1[uindex1];
        Tneta = nt*in2[uindex2];

        double t2 = cos(th);
        double t3 = cos(ps);
        double t4 = cos(om);
        double t5 = cos(ne);
        double t6 = sin(th);
        double t7 = sin(ps);
        double t8 = sin(om);
        double t9 = sin(ne);
        double t10 = mf+mt;
        double t11 = mf*rf;
        double t12 = mt*rt;
        double t13 = rf+rt;
        double t14 = Ifxx*2.0;
        double t15 = Ifxx*4.0;
        double t16 = Ifyy*4.0;
        double t17 = Ifzz*2.0;
        double t18 = Ifzz*4.0;
        double t19 = Itxx*2.0;
        double t20 = Itxx*3.0;
        double t21 = Itxx*4.0;
        double t22 = Ityy*2.0;
        double t23 = Ityy*3.0;
        double t24 = Ityy*4.0;
        double t25 = Itzz*2.0;
        double t26 = Itzz*4.0;
        double t27 = Itzz*Itzz;
        double t28 = Iwxx*2.0;
        double t29 = Iwxx*4.0;
        double t30 = Iwyy*4.0;
        double t31 = Iwzz*2.0;
        double t32 = Iwzz*4.0;
        double t33 = mf*2.0;
        double t34 = mt*2.0;
        double t35 = mt*4.0;
        double t36 = mw*2.0;
        double t37 = mw*4.0;
        double t38 = rf*rf;

        double t40 = rw*rw;
        double t41 = th*2.0;
        double t42 = ps*2.0;
        double t43 = om*2.0;
        double t44 = ne*2.0;
        double t45 = dps*2.0;
        double t46 = dom*2.0;
        double t47 = dph*dph;
        double t48 = dth*dth;
        double t49 = dom*dom;
        double t73 = Ifxx*8.0;
        double t74 = Ifyy*8.0;
        double t75 = -Ifzz;
        double t78 = Ifzz*8.0;
        double t79 = -Ityy;
        double t83 = -Itzz;
        double t86 = Itzz*8.0;
        double t89 = Iwxx*8.0;
        double t90 = Iwyy*8.0;
        double t91 = -Iwzz;
        double t94 = Iwzz*8.0;
        double t97 = -Tneta;
        double t50 = cos(t41);
        double t51 = cos(t42);
        double t52 = cos(t43);
        double t53 = cos(t44);
        double t54 = t2*t2;
        double t55 = t3*t3;
        double t56 = t4*t4;
        double t57 = t4*t4*t4;
        double t58 = t5*t5;
        double t59 = t11*2.0;
        double t60 = t11*4.0;
        double t61 = t12*2.0;
        double t62 = sin(t41);
        double t63 = sin(t42);
        double t64 = sin(t43);
        double t65 = sin(t44);
        double t66 = t6*t6;
        double t67 = t7*t7;
        double t68 = t8*t8;
        double t69 = t9*t9;
        double t70 = mw+t10;
        double t71 = -t14;
        double t72 = -t15;
        double t76 = -t17;
        double t77 = -t18;
        double t80 = -t22;
        double t81 = -t23;
        double t82 = -t24;
        double t84 = -t25;
        double t85 = -t26;
        double t87 = -t28;
        double t88 = -t29;
        double t92 = -t31;
        double t93 = -t32;
        double t95 = rw*t11;
        double t96 = t6*dph;
        double t98 = rf*t10;
        double t99 = mt*t13;
        double t100 = t5*t6;
        double t101 = t40*2.0;
        double t102 = t6*t9;
        double t103 = t11*8.0;
        double t104 = t12*8.0;
        double t105 = rf*t11;
        double t106 = mt*t38;
        double t107 = mf*t40;
        double t108 = rt*t12;
        double t109 = mt*t40;
        double t110 = mw*t40;
        double t111 = t13*t13;
        double t112 = -t78;
        double t113 = -t86;
        double t114 = -t94;
        double t119 = t13*t34;
        double t120 = t13*t35;
        double t121 = -t48;
        double t123 = Itxx+t79;
        double t124 = Iwxx+t91;
        double t142 = t2*t5*t8;
        double t145 = t2*t8*t9;
        double t152 = t34+t36;
        double t153 = t35+t37;
        double t154 = Ifxx*t2*t4*t8;

        double t157 = Iwxx*t2*t3*t7;

        double t180 = t2*t4*t8*t75;
        double t181 = t2*t4*t8*t83;
        double t182 = t2*t3*t7*t91;
        double t188 = t10*2.0+t36;
        double t115 = rf*t61;
        double t116 = t96*2.0;
        double t117 = t96*4.0;
        double t118 = t98*2.0;
        double t122 = t46*t96;
        double t125 = rw*t70;
        double t126 = rw*t99;
        double t127 = t96+dps;
        double t128 = t96+dom;
        double t129 = t50*3.0;
        double t130 = rf*t59;
        double t131 = t105*3.0;
        double t132 = rf*t60;
        double t133 = -t96;
        double t134 = Itxx*t58;
        double t135 = Ityy*t58;
        double t136 = Iwxx*t55;
        double t137 = Ifzz*t68;
        double t138 = Itxx*t69;
        double t139 = Ityy*t69;
        double t140 = Itzz*t68;
        double t141 = Iwzz*t67;
        double t143 = t98*8.0;
        double t144 = t99*8.0;
        double t146 = t13*t99;
        double t149 = rf*t103;
        double t150 = t16*t66;
        double t151 = t30*t66;
        double t159 = t34*t111;
        double t161 = t35*t111;
        double t162 = t47*t50;
        double t163 = t47*t54;
        double t165 = t19+t80;
        double t166 = t20+t81;
        double t167 = t21+t82;
        double t168 = t123*t123;
        double t169 = t28+t92;
        double t170 = t29+t93;
        double t171 = t40*t70;
        double t172 = t45+t96;
        double t173 = t46+t96;
        double t174 = -t142;
        double t175 = t38+t101;
        double t177 = t12+t98;
        double t178 = t11+t99;
        double t187 = t89+t114;
        double t189 = t52*t54*2.0;
        double t191 = t18*t54*t56;
        double t192 = t26*t54*t56;
        double t193 = t32*t54*t55;
        double t194 = t40*t152;
        double t195 = t40*t153;
        double t196 = t15*t54*t68;
        double t197 = t29*t54*t67;
        double t200 = t59+t119;
        double t201 = t60+t120;
        double t202 = t53*t123;
        double t203 = t51*t124;
        double t204 = t63*t124;
        double t208 = t8*t53*t62*4.0;
        double t219 = t100+t145;
        double t224 = t6*t65*t123;
        double t230 = t5*t9*t49*t123;
        double t231 = t3*t7*t48*t124;
        double t234 = t40*t188;

        double t254 = t2*t4*t65*t123*2.0;
        double t255 = t62*t65*t123*dph;
        double t258 = t4*t62*t65*t123;
        double t259 = t2*t64*t65*t123;

        double t291 = t5*t9*t56*t121*t123;
        double t295 = t2*t5*t9*t56*t83*t123;
        double t147 = t117+dps;
        double t148 = -t129;
        double t160 = t146*3.0;

        double t176 = t13*t144;
        double t179 = t171*4.0;
        double t183 = t96*t127;
        double t184 = t46+t116+dps;
        double t185 = -t162;
        double t186 = Iwyy+t171;
        double t190 = t178*t178;
        double t198 = mf*t175;
        double t199 = t61+t118;
        double t205 = t194*4.0;
        double t206 = t33*t175;
        double t209 = t4*t178;
        double t210 = t202*2.0;
        double t211 = t202*4.0;
        double t212 = t203*4.0;
        double t213 = t53*t166;
        double t214 = rw*t4*t177;
        double t215 = t63*t169;
        double t217 = t105+t146;
        double t218 = Itzz+t202;
        double t220 = t104+t143;
        double t221 = t103+t144;

        double t225 = -t202;
        double t229 = t51*t169*2.0;
        double t233 = t224*4.0;

        double t236 = g*t2*t8*t178;
        double t240 = rw*t4*t200;
        double t241 = rw*t8*t200;
        double t242 = rw*t8*t201;
        double t243 = t130+t159;
        double t244 = rw*t201*dth*dom;
        double t246 = t102+t174;
        double t248 = t4*t224;
        double t250 = t219*t219;
        double t251 = rw*t8*t49*t177;
        double t253 = t54*t204*dps;
        double t257 = t2*t4*t65*t165;
        double t260 = t65*t96*t165*dth;
        double t262 = t9*t123*t174;

        double t264 = t63*t187*dth*dps;
        double t265 = t65*t123*t128;
        double t269 = t56*t65*t165*dth;
        double t270 = t259*2.0;
        double t275 = t3*t7*t124*t163;
        double t276 = t259*dph;
        double t280 = Ityy*t4*t9*t219;
        double t282 = rw*t8*t133*t177*dom;
        double t297 = t8*t50*t65*t167*dph;
        double t306 = t56*t58*t69*t168;
        double t316 = t2*t51*t170*t172*dph;
        double t327 = Itzz*rw*t2*t5*t9*t57*t123*t178;
        double t335 = Ifxx+t105+t106+t108+t115+t134+t139;
        double t336 = Ifyy+t105+t106+t108+t115+t135+t138;
        double t207 = t198*4.0;
        double t216 = rw*t209;
        double t223 = t30+t179;
        double t226 = -t210;

        double t228 = t213*2.0;
        double t232 = t48+t183;
        double t238 = -t213;
        double t245 = Itzz+t225;

        double t252 = t4*t217;
        double t256 = t25+t210;
        double t261 = -t236;

        double t268 = t246*t246;
        double t271 = -t251;
        double t272 = rw*t6*t8*t199*2.0;
        double t273 = fcoeff+t253;
        double t277 = t24*t250;
        double t278 = rw*t2*t220*dth;

        double t284 = t40*t56*t190;
        double t285 = t66*t242;
        double t286 = t2*t8*t218*4.0;
        double t287 = t125+t209;
        double t288 = t148+t189+1.0;
        double t289 = t6*t45*t241;
        double t290 = -t269;
        double t292 = -t275;
        double t296 = t49+t122+t185;
        double t305 = rw*t4*t66*t221;
        double t307 = Itxx*t4*t5*t246;
        double t308 = -t297;

        double t314 = t52*t54*t243;
        double t317 = rw*t2*t4*t147*t201;

        double t328 = rw*t2*t4*t184*t221*dph;
        double t329 = t146+t194+t198;
        double t331 = t186*t295;
        double t337 = t186+t203+t214;
        double t349 = t47*t123*t219*t246;
        double t353 = t56*t335;
        double t354 = Itxx+Ityy+t14+t76+t84+t202+t243;
        double t360 = t75+t83+t335;
        double t365 = t15+t19+t22+t77+t85+t132+t161+t210;
        double t369 = t186*t336;
        double t376 = t21+t24+t73+t112+t113+t149+t176+t211;
        double t239 = -t228;
        double t247 = t223*dps;
        double t267 = t25+t226;
        double t283 = t4*t245*dth;
        double t293 = t21*t268;
        double t294 = t2*t8*t245;
        double t298 = t4*t6*t256;
        double t301 = t273*dph*4.0;

        double t303 = rw*t287;
        double t313 = -t307;
        double t315 = t65*t288;
        double t318 = g*t6*t287*8.0;
        double t320 = t4*t128*t256;
        double t321 = -t314;

        double t323 = rw*t200*t232;
        double t324 = t95+t126+t252;
        double t326 = -dps*(t212-t223);
        double t332 = t65*t123*t296;

        double t339 = t254+t272;
        double t340 = t50*t329;
        double t341 = t255+t278;
        double t343 = t16+t19+t22+t132+t161+t226;

        double t355 = t2*t337*dph*dth;
        double t364 = t64*t354;
        double t371 = t64*t360;
        double t372 = t216+t336;
        double t381 = (t2*t52*t365*dph*dth)/4.0;
        double t382 = t64*t376*dth*dom;

        double t385 = t2*t52*t173*t365*dph;
        double t387 = t2*t52*t365*(t96-dom);
        double t403 = Itxx+Ityy+t16+t30+t71+t76+t84+t87+t92+t159+t195+t206+t238;

        double t415 = t107+t109+t110+t136+t137+t140+t141+t240+t353;
        double t299 = t294*dph;
        double t300 = -t283;
        double t304 = t4*t267*dom;
        double t309 = t4*t6*t267*2.0;
        double t310 = Iwyy+t303;
        double t333 = t2*t8*t324;
        double t334 = t224+t294;
        double t344 = -t340;

        double t350 = t343*dom;
        double t351 = t259+t298;
        double t352 = t49*t339;
        double t361 = (t4*t341*dph)/4.0;
        double t367 = (t2*t343*dph*dth)/4.0;
        double t368 = t260+t323;
        double t370 = -t123*dph*(t208-t315);
        double t373 = t244+t332;
        double t374 = t54*t364;

        double t377 = t6*t372;

        double t390 = (t364*(t48-t163))/4.0;
        double t392 = t276+t290+t320;

        double t401 = t215+t242+t364;
        double t402 = t204+t241+t371;
        double t408 = t96*t403;
        double t419 = t19+t22+t72+t74+t77+t85+t88+t90+t93+t161+t205+t207+t229+t239;
        double t423 = Itzz*t2*t4*t415;
        double t431 = t231+t271+t282+t292+t355+Tneta;
        double t433 = t284*t415;

        double t440 = t336*t415;

        double t319 = t310*t310;
        double t342 = t334*dph*dom;
        double t346 = t4*t9*t100*t123*t310;
        double t347 = t6*t216*t310;
        double t357 = t351*dph;
        double t358 = -t352;

        double t378 = -t374;
        double t379 = (t8*t368)/2.0;
        double t380 = t8*t373*4.0;
        double t384 = t265+t299+t300;
        double t386 = t6*t310*t336;
        double t395 = t392*dne*4.0;
        double t397 = t262+t377;
        double t404 = t6*t401;
        double t405 = t2*t402;
        double t424 = t96*t419;
        double t437 = t154+t157+t180+t181+t182+t280+t313+t333;

        double t448 = t247+t350+t408;
        double t449 = dne*(t370+dom*(t233-t286)+dth*(t270-t309));

        double t481 = t131+t150+t151+t160+t191+t192+t193+t196+t197+t234+t277+t293+t305+t321+t344;

        double t330 = t66*t319;
        double t356 = -t347;

        double t363 = t8*t83*t347;
        double t388 = t384*dne;

        double t394 = t8*t83*t386;
        double t398 = -t395;
        double t399 = dth*(t304-t357)*(-1.0/2.0);
        double t400 = Itzz*t8*t397;
        double t407 = t186*t397;
        double t409 = t4*t5*t9*t123*t397;
        double t410 = t216*t397;
        double t412 = t258+t285+t378;
        double t416 = t257+t404;
        double t417 = t6*t310*t397;

        double t426 = Itzz*t2*t4*(t248-t405)*(-1.0/2.0);
        double t429 = t4*t5*t9*t123*(t248-t405)*(-1.0/2.0);
        double t430 = t216*(t248-t405)*(-1.0/2.0);
        double t439 = Itzz*t8*t437;

        double t442 = t347*t415;
        double t444 = t336*(t248-t405)*(-1.0/2.0);
        double t446 = t4*t5*t9*t123*t437;
        double t447 = t216*t437;
        double t450 = t2*t448*dph*2.0;
        double t457 = t336*t437;
        double t458 = t397*t415;

        double t463 = t326+t350+t424;
        double t475 = t397*t437;

        double t482 = (t437*(t248-t405))/2.0;
        double t485 = (Itzz*t8*t481)/4.0;
        double t487 = (t186*t481)/4.0;
        double t489 = (t4*t5*t9*t123*t481)/4.0;
        double t490 = (t216*t481)/4.0;

        double t512 = (t415*t481)/4.0;
        double t411 = -t410;
        double t413 = t412*dph*2.0;
        double t414 = t295+t400;
        double t421 = t48*t416;

        double t438 = -Itzz*t8*(t347-t407);
        double t445 = t356*t415;
        double t454 = t346+t430;
        double t464 = t2*t463;
        double t466 = t346+t447;
        double t469 = t230+t291+t342+t349+t399+Tomega;

        double t478 = t409+t444;
        double t479 = ((t248-t405)*(t347-t407))/2.0;
        double t484 = t423+t439;
        double t488 = -t487;
        double t491 = -t490;
        double t492 = t409+t457;
        double t498 = t429+t458;
        double t508 = t446+t458;
        double t513 = t426+t485;
        double t516 = t97+t261+t361+t367+t379+t381+t388+t390;

        double t525 = t4*t5*t9*t123*(t489-(t397*(t248-t405))/2.0);
        double t527 = t475+t489;
        double t537 = t264+t316+t318+t328+t380+t382+t385+t398+t450;
        double t539 = t482+t512;

        double t422 = t8*t83*t414;
        double t425 = t186*t414;
        double t427 = t216*t414;
        double t428 = t289+t413;
        double t451 = t386+t411;
        double t456 = t331+t438;
        double t459 = t216*t454;
        double t470 = t4*t5*t9*t123*t466;
        double t471 = t216*t466;
        double t476 = t6*t310*t466;
        double t483 = t186*t478;
        double t486 = Itzz*t2*t4*t484;
        double t495 = t186*t492;
        double t496 = t6*t310*t492;
        double t497 = t330+t488;
        double t499 = t308+t317+t387+t464;
        double t504 = t186*t498;
        double t506 = t216*t498;
        double t509 = t6*t310*t498;
        double t511 = t186*t508;
        double t514 = t417+t491;
        double t515 = Itzz*t8*t513;
        double t520 = t397*t508;
        double t529 = t186*t527;
        double t531 = t216*t527;

        double t541 = t83*t539;
        double t542 = t216*t539;
        double t544 = t336*t539;
        double t432 = t428*dom;
        double t452 = t363+t425;
        double t453 = Itzz*t8*t451;

        double t465 = t394+t427;

        double t493 = t415*t451;

        double t500 = t4*t5*t9*t123*t497;
        double t502 = t499*dth;

        double t510 = t336*t497;
        double t517 = t4*t5*t9*t123*t514;
        double t518 = t216*t514;

        double t533 = -t531;
        double t534 = t445+t511;
        double t543 = -t542;
        double t545 = -t544;

        double t559 = t486+t515+t541;
        double t434 = -t432;
        double t467 = t327+t453;
        double t503 = -t500;

        double t549 = t496+t533;

        double t556 = t509+t543;
        double t560 = t216*t559;

        double t567 = Itzz*t8*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))*8.0;
        double t568 = t520+t525+t545;
        double t538 = t479+t503;
        double t548 = t301+t358+t421+t434+t449+t502;

        double t557 = t216*t556;

        double t569 = t186*t568;

        double t576 = -1.0/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double t577 = -1.0/(t26*(t557-t569+t6*t310*(t470-t493))+t8*t85*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+t2*t4*t26*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double t579 = -1.0/(-t567+t86*(t557-t569+t6*t310*(t470-t493))+t2*t4*t86*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));

        double et1 = (t548*(Itzz*(t433+t186*(t306-t440))-t27*t68*(t284-t369)))/(t26*(t557-t569+t6*t310*(t470-t493))+t8*t85*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+t2*t4*t26*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t431*t576*(Itzz*(t470-t493)-Itzz*t8*t465);
        double et2 = (t469*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t537*t579*(Itzz*(t471-t495)+t2*t4*t8*t27*(t284-t369));
        double et3 = (t516*(Itzz*t534+t8*t83*t452))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et4 = (t537*(Itzz*(-t510+t518+t397*(t347-t407))+t27*t54*t56*(t284-t369)))/(-t567+t86*(t557-t569+t6*t310*(t470-t493))+t2*t4*t86*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+(t516*(Itzz*t538+Itzz*t2*t4*t452))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et5 = (t431*(Itzz*(t517-(t451*(t248-t405))/2.0)+t2*t4*t83*t465))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t469*t576*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483));
        double et6 = (t548*(Itzz*(t459-t483)+t2*t4*t8*t27*(t284-t369)))/(t26*(t557-t569+t6*t310*(t470-t493))+t8*t85*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+t2*t4*t26*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et7 = (t431*(Itzz*t568+t8*t83*(t397*t414-t336*t513)+Itzz*t2*t4*(t336*t484+t4*t5*t9*t123*t414)))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t516*t576*(t560+t6*t310*(t422+Itzz*t508))+t537*t579*(Itzz*t549+Itzz*t2*t4*t467)+t548*t577*(Itzz*(t506+t6*t310*(t306-t440))-t8*t83*t467);
        double et8 = (t469*(Itzz*t8*t549-Itzz*t2*t4*(t506+t6*t310*(t306-t440))))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et9 = (t516*(t186*t559+t330*(Itzz*t107+Itzz*t109+Itzz*t110+Itzz*t136+Itzz*t137+Itzz*t141+Ifxx*Itzz*t56+Itzz*t56*t105+Itzz*t56*t106+Itzz*t56*t108+Itzz*t56*t134+Itzz*t56*t139+t4*t25*t95+rf*t12*t25*t56+rw*t4*t12*t25+mt*rf*rw*t4*t25)))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et10 = (t537*(Itzz*(t476-t529)-Itzz*t2*t4*t456))/(-t567+t86*(t557-t569+t6*t310*(t470-t493))+t2*t4*t86*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t431*t576*(t560+t6*t310*(t422+Itzz*t498))+t469*t576*(Itzz*t8*(t476-t529)+Itzz*t2*t4*(t442-t504))+t548*t577*(Itzz*(t442-t504)-t8*t83*t456);
        double et11 = t516*t576*(Itzz*t8*t538+Itzz*t2*t4*t534)+(t469*(t557-t569+t6*t310*(t470-t493)))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t537*t579*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t471-t495))+t548*t577*(Itzz*t8*(t459-t483)+Itzz*t2*t4*(t433+t186*(t306-t440)));
        double et12 = t431*t576*(Itzz*t8*(t517-(t451*(t248-t405))/2.0)-Itzz*t2*t4*(t470-t493));

        d_dTH[index] = et4+et5+et6;
        
        index = index + num_threads;
        index1 = (grid_size[0] > 1) ? index : 0;
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

    int index1 = (grid_size[0] > 1) ? index : 0;
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
        th = TH[index1];
        om = OM[index2];
        dph = dPH[index3];
        dth = dTH[index4];
        dps = dPS[index5];
        dom = dOM[index6];
        dne = dNE[index7];

        Tomega = nw*in1[uindex1];
        Tneta = nt*in2[uindex2];

        double t2 = cos(th);
        double t3 = cos(ps);
        double t4 = cos(om);
        double t5 = cos(ne);
        double t6 = sin(th);
        double t7 = sin(ps);
        double t8 = sin(om);
        double t9 = sin(ne);
        double t10 = mf+mt;
        double t11 = mf*rf;
        double t12 = mt*rt;
        double t13 = rf+rt;
        double t14 = Ifxx*2.0;
        double t15 = Ifxx*4.0;
        double t16 = Ifyy*4.0;
        double t17 = Ifzz*2.0;
        double t18 = Ifzz*4.0;
        double t19 = Itxx*2.0;
        double t20 = Itxx*3.0;
        double t21 = Itxx*4.0;
        double t22 = Ityy*2.0;
        double t23 = Ityy*3.0;
        double t24 = Ityy*4.0;
        double t25 = Itzz*2.0;
        double t26 = Itzz*4.0;
        double t27 = Itzz*Itzz;
        double t28 = Iwxx*2.0;
        double t29 = Iwxx*4.0;
        double t30 = Iwyy*4.0;
        double t31 = Iwzz*2.0;
        double t32 = Iwzz*4.0;
        double t33 = mf*2.0;
        double t34 = mt*2.0;
        double t35 = mt*4.0;
        double t36 = mw*2.0;
        double t37 = mw*4.0;
        double t38 = rf*rf;

        double t40 = rw*rw;
        double t41 = th*2.0;
        double t42 = ps*2.0;
        double t43 = om*2.0;
        double t44 = ne*2.0;
        double t45 = dps*2.0;
        double t46 = dom*2.0;
        double t47 = dph*dph;
        double t48 = dth*dth;
        double t49 = dom*dom;
        double t73 = Ifxx*8.0;
        double t74 = Ifyy*8.0;
        double t75 = -Ifzz;
        double t78 = Ifzz*8.0;
        double t79 = -Ityy;
        double t83 = -Itzz;
        double t86 = Itzz*8.0;
        double t89 = Iwxx*8.0;
        double t90 = Iwyy*8.0;
        double t91 = -Iwzz;
        double t94 = Iwzz*8.0;
        double t97 = -Tneta;
        double t50 = cos(t41);
        double t51 = cos(t42);
        double t52 = cos(t43);
        double t53 = cos(t44);
        double t54 = t2*t2;
        double t55 = t3*t3;
        double t56 = t4*t4;
        double t57 = t4*t4*t4;
        double t58 = t5*t5;
        double t59 = t11*2.0;
        double t60 = t11*4.0;
        double t61 = t12*2.0;
        double t62 = sin(t41);
        double t63 = sin(t42);
        double t64 = sin(t43);
        double t65 = sin(t44);
        double t66 = t6*t6;
        double t67 = t7*t7;
        double t68 = t8*t8;
        double t69 = t9*t9;
        double t70 = mw+t10;
        double t71 = -t14;
        double t72 = -t15;
        double t76 = -t17;
        double t77 = -t18;
        double t80 = -t22;
        double t81 = -t23;
        double t82 = -t24;
        double t84 = -t25;
        double t85 = -t26;
        double t87 = -t28;
        double t88 = -t29;
        double t92 = -t31;
        double t93 = -t32;
        double t95 = rw*t11;
        double t96 = t6*dph;
        double t98 = rf*t10;
        double t99 = mt*t13;
        double t100 = t5*t6;
        double t101 = t40*2.0;
        double t102 = t6*t9;
        double t103 = t11*8.0;
        double t104 = t12*8.0;
        double t105 = rf*t11;
        double t106 = mt*t38;
        double t107 = mf*t40;
        double t108 = rt*t12;
        double t109 = mt*t40;
        double t110 = mw*t40;
        double t111 = t13*t13;
        double t112 = -t78;
        double t113 = -t86;
        double t114 = -t94;
        double t119 = t13*t34;
        double t120 = t13*t35;
        double t121 = -t48;
        double t123 = Itxx+t79;
        double t124 = Iwxx+t91;
        double t142 = t2*t5*t8;
        double t145 = t2*t8*t9;
        double t152 = t34+t36;
        double t153 = t35+t37;
        double t154 = Ifxx*t2*t4*t8;

        double t157 = Iwxx*t2*t3*t7;

        double t180 = t2*t4*t8*t75;
        double t181 = t2*t4*t8*t83;
        double t182 = t2*t3*t7*t91;
        double t188 = t10*2.0+t36;
        double t115 = rf*t61;
        double t116 = t96*2.0;
        double t117 = t96*4.0;
        double t118 = t98*2.0;
        double t122 = t46*t96;
        double t125 = rw*t70;
        double t126 = rw*t99;
        double t127 = t96+dps;
        double t128 = t96+dom;
        double t129 = t50*3.0;
        double t130 = rf*t59;
        double t131 = t105*3.0;
        double t132 = rf*t60;
        double t133 = -t96;
        double t134 = Itxx*t58;
        double t135 = Ityy*t58;
        double t136 = Iwxx*t55;
        double t137 = Ifzz*t68;
        double t138 = Itxx*t69;
        double t139 = Ityy*t69;
        double t140 = Itzz*t68;
        double t141 = Iwzz*t67;
        double t143 = t98*8.0;
        double t144 = t99*8.0;
        double t146 = t13*t99;
        double t149 = rf*t103;
        double t150 = t16*t66;
        double t151 = t30*t66;
        double t159 = t34*t111;
        double t161 = t35*t111;
        double t162 = t47*t50;
        double t163 = t47*t54;
        double t165 = t19+t80;
        double t166 = t20+t81;
        double t167 = t21+t82;
        double t168 = t123*t123;
        double t169 = t28+t92;
        double t170 = t29+t93;
        double t171 = t40*t70;
        double t172 = t45+t96;
        double t173 = t46+t96;
        double t174 = -t142;
        double t175 = t38+t101;
        double t177 = t12+t98;
        double t178 = t11+t99;
        double t187 = t89+t114;
        double t189 = t52*t54*2.0;
        double t191 = t18*t54*t56;
        double t192 = t26*t54*t56;
        double t193 = t32*t54*t55;
        double t194 = t40*t152;
        double t195 = t40*t153;
        double t196 = t15*t54*t68;
        double t197 = t29*t54*t67;
        double t200 = t59+t119;
        double t201 = t60+t120;
        double t202 = t53*t123;
        double t203 = t51*t124;
        double t204 = t63*t124;
        double t208 = t8*t53*t62*4.0;
        double t219 = t100+t145;
        double t224 = t6*t65*t123;
        double t230 = t5*t9*t49*t123;
        double t231 = t3*t7*t48*t124;
        double t234 = t40*t188;

        double t254 = t2*t4*t65*t123*2.0;
        double t255 = t62*t65*t123*dph;
        double t258 = t4*t62*t65*t123;
        double t259 = t2*t64*t65*t123;

        double t291 = t5*t9*t56*t121*t123;
        double t295 = t2*t5*t9*t56*t83*t123;
        double t147 = t117+dps;
        double t148 = -t129;
        double t160 = t146*3.0;

        double t176 = t13*t144;
        double t179 = t171*4.0;
        double t183 = t96*t127;
        double t184 = t46+t116+dps;
        double t185 = -t162;
        double t186 = Iwyy+t171;
        double t190 = t178*t178;
        double t198 = mf*t175;
        double t199 = t61+t118;
        double t205 = t194*4.0;
        double t206 = t33*t175;
        double t209 = t4*t178;
        double t210 = t202*2.0;
        double t211 = t202*4.0;
        double t212 = t203*4.0;
        double t213 = t53*t166;
        double t214 = rw*t4*t177;
        double t215 = t63*t169;
        double t217 = t105+t146;
        double t218 = Itzz+t202;
        double t220 = t104+t143;
        double t221 = t103+t144;

        double t225 = -t202;
        double t229 = t51*t169*2.0;
        double t233 = t224*4.0;

        double t236 = g*t2*t8*t178;
        double t240 = rw*t4*t200;
        double t241 = rw*t8*t200;
        double t242 = rw*t8*t201;
        double t243 = t130+t159;
        double t244 = rw*t201*dth*dom;
        double t246 = t102+t174;
        double t248 = t4*t224;
        double t250 = t219*t219;
        double t251 = rw*t8*t49*t177;
        double t253 = t54*t204*dps;
        double t257 = t2*t4*t65*t165;
        double t260 = t65*t96*t165*dth;
        double t262 = t9*t123*t174;

        double t264 = t63*t187*dth*dps;
        double t265 = t65*t123*t128;
        double t269 = t56*t65*t165*dth;
        double t270 = t259*2.0;
        double t275 = t3*t7*t124*t163;
        double t276 = t259*dph;
        double t280 = Ityy*t4*t9*t219;
        double t282 = rw*t8*t133*t177*dom;
        double t297 = t8*t50*t65*t167*dph;
        double t306 = t56*t58*t69*t168;
        double t316 = t2*t51*t170*t172*dph;
        double t327 = Itzz*rw*t2*t5*t9*t57*t123*t178;
        double t335 = Ifxx+t105+t106+t108+t115+t134+t139;
        double t336 = Ifyy+t105+t106+t108+t115+t135+t138;
        double t207 = t198*4.0;
        double t216 = rw*t209;
        double t223 = t30+t179;
        double t226 = -t210;

        double t228 = t213*2.0;
        double t232 = t48+t183;
        double t238 = -t213;
        double t245 = Itzz+t225;

        double t252 = t4*t217;
        double t256 = t25+t210;
        double t261 = -t236;

        double t268 = t246*t246;
        double t271 = -t251;
        double t272 = rw*t6*t8*t199*2.0;
        double t273 = fcoeff+t253;
        double t277 = t24*t250;
        double t278 = rw*t2*t220*dth;

        double t284 = t40*t56*t190;
        double t285 = t66*t242;
        double t286 = t2*t8*t218*4.0;
        double t287 = t125+t209;
        double t288 = t148+t189+1.0;
        double t289 = t6*t45*t241;
        double t290 = -t269;
        double t292 = -t275;
        double t296 = t49+t122+t185;
        double t305 = rw*t4*t66*t221;
        double t307 = Itxx*t4*t5*t246;
        double t308 = -t297;

        double t314 = t52*t54*t243;
        double t317 = rw*t2*t4*t147*t201;

        double t328 = rw*t2*t4*t184*t221*dph;
        double t329 = t146+t194+t198;
        double t331 = t186*t295;
        double t337 = t186+t203+t214;
        double t349 = t47*t123*t219*t246;
        double t353 = t56*t335;
        double t354 = Itxx+Ityy+t14+t76+t84+t202+t243;
        double t360 = t75+t83+t335;
        double t365 = t15+t19+t22+t77+t85+t132+t161+t210;
        double t369 = t186*t336;
        double t376 = t21+t24+t73+t112+t113+t149+t176+t211;
        double t239 = -t228;
        double t247 = t223*dps;
        double t267 = t25+t226;
        double t283 = t4*t245*dth;
        double t293 = t21*t268;
        double t294 = t2*t8*t245;
        double t298 = t4*t6*t256;
        double t301 = t273*dph*4.0;

        double t303 = rw*t287;
        double t313 = -t307;
        double t315 = t65*t288;
        double t318 = g*t6*t287*8.0;
        double t320 = t4*t128*t256;
        double t321 = -t314;

        double t323 = rw*t200*t232;
        double t324 = t95+t126+t252;
        double t326 = -dps*(t212-t223);
        double t332 = t65*t123*t296;

        double t339 = t254+t272;
        double t340 = t50*t329;
        double t341 = t255+t278;
        double t343 = t16+t19+t22+t132+t161+t226;

        double t355 = t2*t337*dph*dth;
        double t364 = t64*t354;
        double t371 = t64*t360;
        double t372 = t216+t336;
        double t381 = (t2*t52*t365*dph*dth)/4.0;
        double t382 = t64*t376*dth*dom;

        double t385 = t2*t52*t173*t365*dph;
        double t387 = t2*t52*t365*(t96-dom);
        double t403 = Itxx+Ityy+t16+t30+t71+t76+t84+t87+t92+t159+t195+t206+t238;

        double t415 = t107+t109+t110+t136+t137+t140+t141+t240+t353;
        double t299 = t294*dph;
        double t300 = -t283;
        double t304 = t4*t267*dom;
        double t309 = t4*t6*t267*2.0;
        double t310 = Iwyy+t303;
        double t333 = t2*t8*t324;
        double t334 = t224+t294;
        double t344 = -t340;

        double t350 = t343*dom;
        double t351 = t259+t298;
        double t352 = t49*t339;
        double t361 = (t4*t341*dph)/4.0;
        double t367 = (t2*t343*dph*dth)/4.0;
        double t368 = t260+t323;
        double t370 = -t123*dph*(t208-t315);
        double t373 = t244+t332;
        double t374 = t54*t364;

        double t377 = t6*t372;

        double t390 = (t364*(t48-t163))/4.0;
        double t392 = t276+t290+t320;

        double t401 = t215+t242+t364;
        double t402 = t204+t241+t371;
        double t408 = t96*t403;
        double t419 = t19+t22+t72+t74+t77+t85+t88+t90+t93+t161+t205+t207+t229+t239;
        double t423 = Itzz*t2*t4*t415;
        double t431 = t231+t271+t282+t292+t355+Tneta;
        double t433 = t284*t415;

        double t440 = t336*t415;

        double t319 = t310*t310;
        double t342 = t334*dph*dom;
        double t346 = t4*t9*t100*t123*t310;
        double t347 = t6*t216*t310;
        double t357 = t351*dph;
        double t358 = -t352;

        double t378 = -t374;
        double t379 = (t8*t368)/2.0;
        double t380 = t8*t373*4.0;
        double t384 = t265+t299+t300;
        double t386 = t6*t310*t336;
        double t395 = t392*dne*4.0;
        double t397 = t262+t377;
        double t404 = t6*t401;
        double t405 = t2*t402;
        double t424 = t96*t419;
        double t437 = t154+t157+t180+t181+t182+t280+t313+t333;

        double t448 = t247+t350+t408;
        double t449 = dne*(t370+dom*(t233-t286)+dth*(t270-t309));

        double t481 = t131+t150+t151+t160+t191+t192+t193+t196+t197+t234+t277+t293+t305+t321+t344;

        double t330 = t66*t319;
        double t356 = -t347;

        double t363 = t8*t83*t347;
        double t388 = t384*dne;

        double t394 = t8*t83*t386;
        double t398 = -t395;
        double t399 = dth*(t304-t357)*(-1.0/2.0);
        double t400 = Itzz*t8*t397;
        double t407 = t186*t397;
        double t409 = t4*t5*t9*t123*t397;
        double t410 = t216*t397;
        double t412 = t258+t285+t378;
        double t416 = t257+t404;
        double t417 = t6*t310*t397;

        double t426 = Itzz*t2*t4*(t248-t405)*(-1.0/2.0);
        double t429 = t4*t5*t9*t123*(t248-t405)*(-1.0/2.0);
        double t430 = t216*(t248-t405)*(-1.0/2.0);
        double t439 = Itzz*t8*t437;

        double t442 = t347*t415;
        double t444 = t336*(t248-t405)*(-1.0/2.0);
        double t446 = t4*t5*t9*t123*t437;
        double t447 = t216*t437;
        double t450 = t2*t448*dph*2.0;
        double t457 = t336*t437;
        double t458 = t397*t415;

        double t463 = t326+t350+t424;
        double t475 = t397*t437;

        double t482 = (t437*(t248-t405))/2.0;
        double t485 = (Itzz*t8*t481)/4.0;
        double t487 = (t186*t481)/4.0;
        double t489 = (t4*t5*t9*t123*t481)/4.0;
        double t490 = (t216*t481)/4.0;

        double t512 = (t415*t481)/4.0;
        double t411 = -t410;
        double t413 = t412*dph*2.0;
        double t414 = t295+t400;
        double t421 = t48*t416;

        double t438 = -Itzz*t8*(t347-t407);
        double t445 = t356*t415;
        double t454 = t346+t430;
        double t464 = t2*t463;
        double t466 = t346+t447;
        double t469 = t230+t291+t342+t349+t399+Tomega;

        double t478 = t409+t444;
        double t479 = ((t248-t405)*(t347-t407))/2.0;
        double t484 = t423+t439;
        double t488 = -t487;
        double t491 = -t490;
        double t492 = t409+t457;
        double t498 = t429+t458;
        double t508 = t446+t458;
        double t513 = t426+t485;
        double t516 = t97+t261+t361+t367+t379+t381+t388+t390;

        double t525 = t4*t5*t9*t123*(t489-(t397*(t248-t405))/2.0);
        double t527 = t475+t489;
        double t537 = t264+t316+t318+t328+t380+t382+t385+t398+t450;
        double t539 = t482+t512;

        double t422 = t8*t83*t414;
        double t425 = t186*t414;
        double t427 = t216*t414;
        double t428 = t289+t413;
        double t451 = t386+t411;
        double t456 = t331+t438;
        double t459 = t216*t454;
        double t470 = t4*t5*t9*t123*t466;
        double t471 = t216*t466;
        double t476 = t6*t310*t466;
        double t483 = t186*t478;
        double t486 = Itzz*t2*t4*t484;
        double t495 = t186*t492;
        double t496 = t6*t310*t492;
        double t497 = t330+t488;
        double t499 = t308+t317+t387+t464;
        double t504 = t186*t498;
        double t506 = t216*t498;
        double t509 = t6*t310*t498;
        double t511 = t186*t508;
        double t514 = t417+t491;
        double t515 = Itzz*t8*t513;
        double t520 = t397*t508;
        double t529 = t186*t527;
        double t531 = t216*t527;

        double t541 = t83*t539;
        double t542 = t216*t539;
        double t544 = t336*t539;
        double t432 = t428*dom;
        double t452 = t363+t425;
        double t453 = Itzz*t8*t451;

        double t465 = t394+t427;

        double t493 = t415*t451;

        double t500 = t4*t5*t9*t123*t497;
        double t502 = t499*dth;

        double t510 = t336*t497;
        double t517 = t4*t5*t9*t123*t514;
        double t518 = t216*t514;

        double t533 = -t531;
        double t534 = t445+t511;
        double t543 = -t542;
        double t545 = -t544;

        double t559 = t486+t515+t541;
        double t434 = -t432;
        double t467 = t327+t453;
        double t503 = -t500;

        double t549 = t496+t533;

        double t556 = t509+t543;
        double t560 = t216*t559;

        double t567 = Itzz*t8*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))*8.0;
        double t568 = t520+t525+t545;
        double t538 = t479+t503;
        double t548 = t301+t358+t421+t434+t449+t502;

        double t557 = t216*t556;

        double t569 = t186*t568;

        double t576 = -1.0/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double t577 = -1.0/(t26*(t557-t569+t6*t310*(t470-t493))+t8*t85*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+t2*t4*t26*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double t579 = -1.0/(-t567+t86*(t557-t569+t6*t310*(t470-t493))+t2*t4*t86*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));

        double et1 = (t548*(Itzz*(t433+t186*(t306-t440))-t27*t68*(t284-t369)))/(t26*(t557-t569+t6*t310*(t470-t493))+t8*t85*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+t2*t4*t26*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t431*t576*(Itzz*(t470-t493)-Itzz*t8*t465);
        double et2 = (t469*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t537*t579*(Itzz*(t471-t495)+t2*t4*t8*t27*(t284-t369));
        double et3 = (t516*(Itzz*t534+t8*t83*t452))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et4 = (t537*(Itzz*(-t510+t518+t397*(t347-t407))+t27*t54*t56*(t284-t369)))/(-t567+t86*(t557-t569+t6*t310*(t470-t493))+t2*t4*t86*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+(t516*(Itzz*t538+Itzz*t2*t4*t452))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et5 = (t431*(Itzz*(t517-(t451*(t248-t405))/2.0)+t2*t4*t83*t465))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t469*t576*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483));
        double et6 = (t548*(Itzz*(t459-t483)+t2*t4*t8*t27*(t284-t369)))/(t26*(t557-t569+t6*t310*(t470-t493))+t8*t85*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+t2*t4*t26*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et7 = (t431*(Itzz*t568+t8*t83*(t397*t414-t336*t513)+Itzz*t2*t4*(t336*t484+t4*t5*t9*t123*t414)))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t516*t576*(t560+t6*t310*(t422+Itzz*t508))+t537*t579*(Itzz*t549+Itzz*t2*t4*t467)+t548*t577*(Itzz*(t506+t6*t310*(t306-t440))-t8*t83*t467);
        double et8 = (t469*(Itzz*t8*t549-Itzz*t2*t4*(t506+t6*t310*(t306-t440))))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et9 = (t516*(t186*t559+t330*(Itzz*t107+Itzz*t109+Itzz*t110+Itzz*t136+Itzz*t137+Itzz*t141+Ifxx*Itzz*t56+Itzz*t56*t105+Itzz*t56*t106+Itzz*t56*t108+Itzz*t56*t134+Itzz*t56*t139+t4*t25*t95+rf*t12*t25*t56+rw*t4*t12*t25+mt*rf*rw*t4*t25)))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et10 = (t537*(Itzz*(t476-t529)-Itzz*t2*t4*t456))/(-t567+t86*(t557-t569+t6*t310*(t470-t493))+t2*t4*t86*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t431*t576*(t560+t6*t310*(t422+Itzz*t498))+t469*t576*(Itzz*t8*(t476-t529)+Itzz*t2*t4*(t442-t504))+t548*t577*(Itzz*(t442-t504)-t8*t83*t456);
        double et11 = t516*t576*(Itzz*t8*t538+Itzz*t2*t4*t534)+(t469*(t557-t569+t6*t310*(t470-t493)))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t537*t579*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t471-t495))+t548*t577*(Itzz*t8*(t459-t483)+Itzz*t2*t4*(t433+t186*(t306-t440)));
        double et12 = t431*t576*(Itzz*t8*(t517-(t451*(t248-t405))/2.0)-Itzz*t2*t4*(t470-t493));

        d_dPS[index] = et7+et8;
        
        index = index + num_threads;
        index1 = (grid_size[0] > 1) ? index : 0;
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
    
    int index1 = (grid_size[0] > 1) ? index : 0;
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
        th = TH[index1];
        om = OM[index2];
        dph = dPH[index3];
        dth = dTH[index4];
        dps = dPS[index5];
        dom = dOM[index6];
        dne = dNE[index7];

        Tomega = nw*in1[uindex1];
        Tneta = nt*in2[uindex2];

        double t2 = cos(th);
        double t3 = cos(ps);
        double t4 = cos(om);
        double t5 = cos(ne);
        double t6 = sin(th);
        double t7 = sin(ps);
        double t8 = sin(om);
        double t9 = sin(ne);
        double t10 = mf+mt;
        double t11 = mf*rf;
        double t12 = mt*rt;
        double t13 = rf+rt;
        double t14 = Ifxx*2.0;
        double t15 = Ifxx*4.0;
        double t16 = Ifyy*4.0;
        double t17 = Ifzz*2.0;
        double t18 = Ifzz*4.0;
        double t19 = Itxx*2.0;
        double t20 = Itxx*3.0;
        double t21 = Itxx*4.0;
        double t22 = Ityy*2.0;
        double t23 = Ityy*3.0;
        double t24 = Ityy*4.0;
        double t25 = Itzz*2.0;
        double t26 = Itzz*4.0;
        double t27 = Itzz*Itzz;
        double t28 = Iwxx*2.0;
        double t29 = Iwxx*4.0;
        double t30 = Iwyy*4.0;
        double t31 = Iwzz*2.0;
        double t32 = Iwzz*4.0;
        double t33 = mf*2.0;
        double t34 = mt*2.0;
        double t35 = mt*4.0;
        double t36 = mw*2.0;
        double t37 = mw*4.0;
        double t38 = rf*rf;

        double t40 = rw*rw;
        double t41 = th*2.0;
        double t42 = ps*2.0;
        double t43 = om*2.0;
        double t44 = ne*2.0;
        double t45 = dps*2.0;
        double t46 = dom*2.0;
        double t47 = dph*dph;
        double t48 = dth*dth;
        double t49 = dom*dom;
        double t73 = Ifxx*8.0;
        double t74 = Ifyy*8.0;
        double t75 = -Ifzz;
        double t78 = Ifzz*8.0;
        double t79 = -Ityy;
        double t83 = -Itzz;
        double t86 = Itzz*8.0;
        double t89 = Iwxx*8.0;
        double t90 = Iwyy*8.0;
        double t91 = -Iwzz;
        double t94 = Iwzz*8.0;
        double t97 = -Tneta;
        double t50 = cos(t41);
        double t51 = cos(t42);
        double t52 = cos(t43);
        double t53 = cos(t44);
        double t54 = t2*t2;
        double t55 = t3*t3;
        double t56 = t4*t4;
        double t57 = t4*t4*t4;
        double t58 = t5*t5;
        double t59 = t11*2.0;
        double t60 = t11*4.0;
        double t61 = t12*2.0;
        double t62 = sin(t41);
        double t63 = sin(t42);
        double t64 = sin(t43);
        double t65 = sin(t44);
        double t66 = t6*t6;
        double t67 = t7*t7;
        double t68 = t8*t8;
        double t69 = t9*t9;
        double t70 = mw+t10;
        double t71 = -t14;
        double t72 = -t15;
        double t76 = -t17;
        double t77 = -t18;
        double t80 = -t22;
        double t81 = -t23;
        double t82 = -t24;
        double t84 = -t25;
        double t85 = -t26;
        double t87 = -t28;
        double t88 = -t29;
        double t92 = -t31;
        double t93 = -t32;
        double t95 = rw*t11;
        double t96 = t6*dph;
        double t98 = rf*t10;
        double t99 = mt*t13;
        double t100 = t5*t6;
        double t101 = t40*2.0;
        double t102 = t6*t9;
        double t103 = t11*8.0;
        double t104 = t12*8.0;
        double t105 = rf*t11;
        double t106 = mt*t38;
        double t107 = mf*t40;
        double t108 = rt*t12;
        double t109 = mt*t40;
        double t110 = mw*t40;
        double t111 = t13*t13;
        double t112 = -t78;
        double t113 = -t86;
        double t114 = -t94;
        double t119 = t13*t34;
        double t120 = t13*t35;
        double t121 = -t48;
        double t123 = Itxx+t79;
        double t124 = Iwxx+t91;
        double t142 = t2*t5*t8;
        double t145 = t2*t8*t9;
        double t152 = t34+t36;
        double t153 = t35+t37;
        double t154 = Ifxx*t2*t4*t8;

        double t157 = Iwxx*t2*t3*t7;

        double t180 = t2*t4*t8*t75;
        double t181 = t2*t4*t8*t83;
        double t182 = t2*t3*t7*t91;
        double t188 = t10*2.0+t36;
        double t115 = rf*t61;
        double t116 = t96*2.0;
        double t117 = t96*4.0;
        double t118 = t98*2.0;
        double t122 = t46*t96;
        double t125 = rw*t70;
        double t126 = rw*t99;
        double t127 = t96+dps;
        double t128 = t96+dom;
        double t129 = t50*3.0;
        double t130 = rf*t59;
        double t131 = t105*3.0;
        double t132 = rf*t60;
        double t133 = -t96;
        double t134 = Itxx*t58;
        double t135 = Ityy*t58;
        double t136 = Iwxx*t55;
        double t137 = Ifzz*t68;
        double t138 = Itxx*t69;
        double t139 = Ityy*t69;
        double t140 = Itzz*t68;
        double t141 = Iwzz*t67;
        double t143 = t98*8.0;
        double t144 = t99*8.0;
        double t146 = t13*t99;
        double t149 = rf*t103;
        double t150 = t16*t66;
        double t151 = t30*t66;
        double t159 = t34*t111;
        double t161 = t35*t111;
        double t162 = t47*t50;
        double t163 = t47*t54;
        double t165 = t19+t80;
        double t166 = t20+t81;
        double t167 = t21+t82;
        double t168 = t123*t123;
        double t169 = t28+t92;
        double t170 = t29+t93;
        double t171 = t40*t70;
        double t172 = t45+t96;
        double t173 = t46+t96;
        double t174 = -t142;
        double t175 = t38+t101;
        double t177 = t12+t98;
        double t178 = t11+t99;
        double t187 = t89+t114;
        double t189 = t52*t54*2.0;
        double t191 = t18*t54*t56;
        double t192 = t26*t54*t56;
        double t193 = t32*t54*t55;
        double t194 = t40*t152;
        double t195 = t40*t153;
        double t196 = t15*t54*t68;
        double t197 = t29*t54*t67;
        double t200 = t59+t119;
        double t201 = t60+t120;
        double t202 = t53*t123;
        double t203 = t51*t124;
        double t204 = t63*t124;
        double t208 = t8*t53*t62*4.0;
        double t219 = t100+t145;
        double t224 = t6*t65*t123;
        double t230 = t5*t9*t49*t123;
        double t231 = t3*t7*t48*t124;
        double t234 = t40*t188;

        double t254 = t2*t4*t65*t123*2.0;
        double t255 = t62*t65*t123*dph;
        double t258 = t4*t62*t65*t123;
        double t259 = t2*t64*t65*t123;

        double t291 = t5*t9*t56*t121*t123;
        double t295 = t2*t5*t9*t56*t83*t123;
        double t147 = t117+dps;
        double t148 = -t129;
        double t160 = t146*3.0;

        double t176 = t13*t144;
        double t179 = t171*4.0;
        double t183 = t96*t127;
        double t184 = t46+t116+dps;
        double t185 = -t162;
        double t186 = Iwyy+t171;
        double t190 = t178*t178;
        double t198 = mf*t175;
        double t199 = t61+t118;
        double t205 = t194*4.0;
        double t206 = t33*t175;
        double t209 = t4*t178;
        double t210 = t202*2.0;
        double t211 = t202*4.0;
        double t212 = t203*4.0;
        double t213 = t53*t166;
        double t214 = rw*t4*t177;
        double t215 = t63*t169;
        double t217 = t105+t146;
        double t218 = Itzz+t202;
        double t220 = t104+t143;
        double t221 = t103+t144;

        double t225 = -t202;
        double t229 = t51*t169*2.0;
        double t233 = t224*4.0;

        double t236 = g*t2*t8*t178;
        double t240 = rw*t4*t200;
        double t241 = rw*t8*t200;
        double t242 = rw*t8*t201;
        double t243 = t130+t159;
        double t244 = rw*t201*dth*dom;
        double t246 = t102+t174;
        double t248 = t4*t224;
        double t250 = t219*t219;
        double t251 = rw*t8*t49*t177;
        double t253 = t54*t204*dps;
        double t257 = t2*t4*t65*t165;
        double t260 = t65*t96*t165*dth;
        double t262 = t9*t123*t174;

        double t264 = t63*t187*dth*dps;
        double t265 = t65*t123*t128;
        double t269 = t56*t65*t165*dth;
        double t270 = t259*2.0;
        double t275 = t3*t7*t124*t163;
        double t276 = t259*dph;
        double t280 = Ityy*t4*t9*t219;
        double t282 = rw*t8*t133*t177*dom;
        double t297 = t8*t50*t65*t167*dph;
        double t306 = t56*t58*t69*t168;
        double t316 = t2*t51*t170*t172*dph;
        double t327 = Itzz*rw*t2*t5*t9*t57*t123*t178;
        double t335 = Ifxx+t105+t106+t108+t115+t134+t139;
        double t336 = Ifyy+t105+t106+t108+t115+t135+t138;
        double t207 = t198*4.0;
        double t216 = rw*t209;
        double t223 = t30+t179;
        double t226 = -t210;

        double t228 = t213*2.0;
        double t232 = t48+t183;
        double t238 = -t213;
        double t245 = Itzz+t225;

        double t252 = t4*t217;
        double t256 = t25+t210;
        double t261 = -t236;

        double t268 = t246*t246;
        double t271 = -t251;
        double t272 = rw*t6*t8*t199*2.0;
        double t273 = fcoeff+t253;
        double t277 = t24*t250;
        double t278 = rw*t2*t220*dth;

        double t284 = t40*t56*t190;
        double t285 = t66*t242;
        double t286 = t2*t8*t218*4.0;
        double t287 = t125+t209;
        double t288 = t148+t189+1.0;
        double t289 = t6*t45*t241;
        double t290 = -t269;
        double t292 = -t275;
        double t296 = t49+t122+t185;
        double t305 = rw*t4*t66*t221;
        double t307 = Itxx*t4*t5*t246;
        double t308 = -t297;

        double t314 = t52*t54*t243;
        double t317 = rw*t2*t4*t147*t201;

        double t328 = rw*t2*t4*t184*t221*dph;
        double t329 = t146+t194+t198;
        double t331 = t186*t295;
        double t337 = t186+t203+t214;
        double t349 = t47*t123*t219*t246;
        double t353 = t56*t335;
        double t354 = Itxx+Ityy+t14+t76+t84+t202+t243;
        double t360 = t75+t83+t335;
        double t365 = t15+t19+t22+t77+t85+t132+t161+t210;
        double t369 = t186*t336;
        double t376 = t21+t24+t73+t112+t113+t149+t176+t211;
        double t239 = -t228;
        double t247 = t223*dps;
        double t267 = t25+t226;
        double t283 = t4*t245*dth;
        double t293 = t21*t268;
        double t294 = t2*t8*t245;
        double t298 = t4*t6*t256;
        double t301 = t273*dph*4.0;

        double t303 = rw*t287;
        double t313 = -t307;
        double t315 = t65*t288;
        double t318 = g*t6*t287*8.0;
        double t320 = t4*t128*t256;
        double t321 = -t314;

        double t323 = rw*t200*t232;
        double t324 = t95+t126+t252;
        double t326 = -dps*(t212-t223);
        double t332 = t65*t123*t296;

        double t339 = t254+t272;
        double t340 = t50*t329;
        double t341 = t255+t278;
        double t343 = t16+t19+t22+t132+t161+t226;

        double t355 = t2*t337*dph*dth;
        double t364 = t64*t354;
        double t371 = t64*t360;
        double t372 = t216+t336;
        double t381 = (t2*t52*t365*dph*dth)/4.0;
        double t382 = t64*t376*dth*dom;

        double t385 = t2*t52*t173*t365*dph;
        double t387 = t2*t52*t365*(t96-dom);
        double t403 = Itxx+Ityy+t16+t30+t71+t76+t84+t87+t92+t159+t195+t206+t238;

        double t415 = t107+t109+t110+t136+t137+t140+t141+t240+t353;
        double t299 = t294*dph;
        double t300 = -t283;
        double t304 = t4*t267*dom;
        double t309 = t4*t6*t267*2.0;
        double t310 = Iwyy+t303;
        double t333 = t2*t8*t324;
        double t334 = t224+t294;
        double t344 = -t340;

        double t350 = t343*dom;
        double t351 = t259+t298;
        double t352 = t49*t339;
        double t361 = (t4*t341*dph)/4.0;
        double t367 = (t2*t343*dph*dth)/4.0;
        double t368 = t260+t323;
        double t370 = -t123*dph*(t208-t315);
        double t373 = t244+t332;
        double t374 = t54*t364;

        double t377 = t6*t372;

        double t390 = (t364*(t48-t163))/4.0;
        double t392 = t276+t290+t320;

        double t401 = t215+t242+t364;
        double t402 = t204+t241+t371;
        double t408 = t96*t403;
        double t419 = t19+t22+t72+t74+t77+t85+t88+t90+t93+t161+t205+t207+t229+t239;
        double t423 = Itzz*t2*t4*t415;
        double t431 = t231+t271+t282+t292+t355+Tneta;
        double t433 = t284*t415;

        double t440 = t336*t415;

        double t319 = t310*t310;
        double t342 = t334*dph*dom;
        double t346 = t4*t9*t100*t123*t310;
        double t347 = t6*t216*t310;
        double t357 = t351*dph;
        double t358 = -t352;

        double t378 = -t374;
        double t379 = (t8*t368)/2.0;
        double t380 = t8*t373*4.0;
        double t384 = t265+t299+t300;
        double t386 = t6*t310*t336;
        double t395 = t392*dne*4.0;
        double t397 = t262+t377;
        double t404 = t6*t401;
        double t405 = t2*t402;
        double t424 = t96*t419;
        double t437 = t154+t157+t180+t181+t182+t280+t313+t333;

        double t448 = t247+t350+t408;
        double t449 = dne*(t370+dom*(t233-t286)+dth*(t270-t309));

        double t481 = t131+t150+t151+t160+t191+t192+t193+t196+t197+t234+t277+t293+t305+t321+t344;

        double t330 = t66*t319;
        double t356 = -t347;

        double t363 = t8*t83*t347;
        double t388 = t384*dne;

        double t394 = t8*t83*t386;
        double t398 = -t395;
        double t399 = dth*(t304-t357)*(-1.0/2.0);
        double t400 = Itzz*t8*t397;
        double t407 = t186*t397;
        double t409 = t4*t5*t9*t123*t397;
        double t410 = t216*t397;
        double t412 = t258+t285+t378;
        double t416 = t257+t404;
        double t417 = t6*t310*t397;

        double t426 = Itzz*t2*t4*(t248-t405)*(-1.0/2.0);
        double t429 = t4*t5*t9*t123*(t248-t405)*(-1.0/2.0);
        double t430 = t216*(t248-t405)*(-1.0/2.0);
        double t439 = Itzz*t8*t437;

        double t442 = t347*t415;
        double t444 = t336*(t248-t405)*(-1.0/2.0);
        double t446 = t4*t5*t9*t123*t437;
        double t447 = t216*t437;
        double t450 = t2*t448*dph*2.0;
        double t457 = t336*t437;
        double t458 = t397*t415;

        double t463 = t326+t350+t424;
        double t475 = t397*t437;

        double t482 = (t437*(t248-t405))/2.0;
        double t485 = (Itzz*t8*t481)/4.0;
        double t487 = (t186*t481)/4.0;
        double t489 = (t4*t5*t9*t123*t481)/4.0;
        double t490 = (t216*t481)/4.0;

        double t512 = (t415*t481)/4.0;
        double t411 = -t410;
        double t413 = t412*dph*2.0;
        double t414 = t295+t400;
        double t421 = t48*t416;

        double t438 = -Itzz*t8*(t347-t407);
        double t445 = t356*t415;
        double t454 = t346+t430;
        double t464 = t2*t463;
        double t466 = t346+t447;
        double t469 = t230+t291+t342+t349+t399+Tomega;

        double t478 = t409+t444;
        double t479 = ((t248-t405)*(t347-t407))/2.0;
        double t484 = t423+t439;
        double t488 = -t487;
        double t491 = -t490;
        double t492 = t409+t457;
        double t498 = t429+t458;
        double t508 = t446+t458;
        double t513 = t426+t485;
        double t516 = t97+t261+t361+t367+t379+t381+t388+t390;

        double t525 = t4*t5*t9*t123*(t489-(t397*(t248-t405))/2.0);
        double t527 = t475+t489;
        double t537 = t264+t316+t318+t328+t380+t382+t385+t398+t450;
        double t539 = t482+t512;

        double t422 = t8*t83*t414;
        double t425 = t186*t414;
        double t427 = t216*t414;
        double t428 = t289+t413;
        double t451 = t386+t411;
        double t456 = t331+t438;
        double t459 = t216*t454;
        double t470 = t4*t5*t9*t123*t466;
        double t471 = t216*t466;
        double t476 = t6*t310*t466;
        double t483 = t186*t478;
        double t486 = Itzz*t2*t4*t484;
        double t495 = t186*t492;
        double t496 = t6*t310*t492;
        double t497 = t330+t488;
        double t499 = t308+t317+t387+t464;
        double t504 = t186*t498;
        double t506 = t216*t498;
        double t509 = t6*t310*t498;
        double t511 = t186*t508;
        double t514 = t417+t491;
        double t515 = Itzz*t8*t513;
        double t520 = t397*t508;
        double t529 = t186*t527;
        double t531 = t216*t527;

        double t541 = t83*t539;
        double t542 = t216*t539;
        double t544 = t336*t539;
        double t432 = t428*dom;
        double t452 = t363+t425;
        double t453 = Itzz*t8*t451;

        double t465 = t394+t427;

        double t493 = t415*t451;

        double t500 = t4*t5*t9*t123*t497;
        double t502 = t499*dth;

        double t510 = t336*t497;
        double t517 = t4*t5*t9*t123*t514;
        double t518 = t216*t514;

        double t533 = -t531;
        double t534 = t445+t511;
        double t543 = -t542;
        double t545 = -t544;

        double t559 = t486+t515+t541;
        double t434 = -t432;
        double t467 = t327+t453;
        double t503 = -t500;

        double t549 = t496+t533;

        double t556 = t509+t543;
        double t560 = t216*t559;

        double t567 = Itzz*t8*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))*8.0;
        double t568 = t520+t525+t545;
        double t538 = t479+t503;
        double t548 = t301+t358+t421+t434+t449+t502;

        double t557 = t216*t556;

        double t569 = t186*t568;

        double t576 = -1.0/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double t577 = -1.0/(t26*(t557-t569+t6*t310*(t470-t493))+t8*t85*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+t2*t4*t26*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double t579 = -1.0/(-t567+t86*(t557-t569+t6*t310*(t470-t493))+t2*t4*t86*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));

        double et1 = (t548*(Itzz*(t433+t186*(t306-t440))-t27*t68*(t284-t369)))/(t26*(t557-t569+t6*t310*(t470-t493))+t8*t85*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+t2*t4*t26*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t431*t576*(Itzz*(t470-t493)-Itzz*t8*t465);
        double et2 = (t469*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t537*t579*(Itzz*(t471-t495)+t2*t4*t8*t27*(t284-t369));
        double et3 = (t516*(Itzz*t534+t8*t83*t452))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et4 = (t537*(Itzz*(-t510+t518+t397*(t347-t407))+t27*t54*t56*(t284-t369)))/(-t567+t86*(t557-t569+t6*t310*(t470-t493))+t2*t4*t86*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+(t516*(Itzz*t538+Itzz*t2*t4*t452))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et5 = (t431*(Itzz*(t517-(t451*(t248-t405))/2.0)+t2*t4*t83*t465))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t469*t576*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483));
        double et6 = (t548*(Itzz*(t459-t483)+t2*t4*t8*t27*(t284-t369)))/(t26*(t557-t569+t6*t310*(t470-t493))+t8*t85*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+t2*t4*t26*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et7 = (t431*(Itzz*t568+t8*t83*(t397*t414-t336*t513)+Itzz*t2*t4*(t336*t484+t4*t5*t9*t123*t414)))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t516*t576*(t560+t6*t310*(t422+Itzz*t508))+t537*t579*(Itzz*t549+Itzz*t2*t4*t467)+t548*t577*(Itzz*(t506+t6*t310*(t306-t440))-t8*t83*t467);
        double et8 = (t469*(Itzz*t8*t549-Itzz*t2*t4*(t506+t6*t310*(t306-t440))))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et9 = (t516*(t186*t559+t330*(Itzz*t107+Itzz*t109+Itzz*t110+Itzz*t136+Itzz*t137+Itzz*t141+Ifxx*Itzz*t56+Itzz*t56*t105+Itzz*t56*t106+Itzz*t56*t108+Itzz*t56*t134+Itzz*t56*t139+t4*t25*t95+rf*t12*t25*t56+rw*t4*t12*t25+mt*rf*rw*t4*t25)))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et10 = (t537*(Itzz*(t476-t529)-Itzz*t2*t4*t456))/(-t567+t86*(t557-t569+t6*t310*(t470-t493))+t2*t4*t86*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t431*t576*(t560+t6*t310*(t422+Itzz*t498))+t469*t576*(Itzz*t8*(t476-t529)+Itzz*t2*t4*(t442-t504))+t548*t577*(Itzz*(t442-t504)-t8*t83*t456);
        double et11 = t516*t576*(Itzz*t8*t538+Itzz*t2*t4*t534)+(t469*(t557-t569+t6*t310*(t470-t493)))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t537*t579*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t471-t495))+t548*t577*(Itzz*t8*(t459-t483)+Itzz*t2*t4*(t433+t186*(t306-t440)));
        double et12 = t431*t576*(Itzz*t8*(t517-(t451*(t248-t405))/2.0)-Itzz*t2*t4*(t470-t493));

        d_dOM[index] = et9+et10;
        
        index = index + num_threads;
        index1 = (grid_size[0] > 1) ? index : 0;
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

    int index1 = (grid_size[0] > 1) ? index : 0;
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
        th = TH[index1];
        om = OM[index2];
        dph = dPH[index3];
        dth = dTH[index4];
        dps = dPS[index5];
        dom = dOM[index6];
        dne = dNE[index7];

        Tomega = nw*in1[uindex1];
        Tneta = nt*in2[uindex2];

        double t2 = cos(th);
        double t3 = cos(ps);
        double t4 = cos(om);
        double t5 = cos(ne);
        double t6 = sin(th);
        double t7 = sin(ps);
        double t8 = sin(om);
        double t9 = sin(ne);
        double t10 = mf+mt;
        double t11 = mf*rf;
        double t12 = mt*rt;
        double t13 = rf+rt;
        double t14 = Ifxx*2.0;
        double t15 = Ifxx*4.0;
        double t16 = Ifyy*4.0;
        double t17 = Ifzz*2.0;
        double t18 = Ifzz*4.0;
        double t19 = Itxx*2.0;
        double t20 = Itxx*3.0;
        double t21 = Itxx*4.0;
        double t22 = Ityy*2.0;
        double t23 = Ityy*3.0;
        double t24 = Ityy*4.0;
        double t25 = Itzz*2.0;
        double t26 = Itzz*4.0;
        double t27 = Itzz*Itzz;
        double t28 = Iwxx*2.0;
        double t29 = Iwxx*4.0;
        double t30 = Iwyy*4.0;
        double t31 = Iwzz*2.0;
        double t32 = Iwzz*4.0;
        double t33 = mf*2.0;
        double t34 = mt*2.0;
        double t35 = mt*4.0;
        double t36 = mw*2.0;
        double t37 = mw*4.0;
        double t38 = rf*rf;

        double t40 = rw*rw;
        double t41 = th*2.0;
        double t42 = ps*2.0;
        double t43 = om*2.0;
        double t44 = ne*2.0;
        double t45 = dps*2.0;
        double t46 = dom*2.0;
        double t47 = dph*dph;
        double t48 = dth*dth;
        double t49 = dom*dom;
        double t73 = Ifxx*8.0;
        double t74 = Ifyy*8.0;
        double t75 = -Ifzz;
        double t78 = Ifzz*8.0;
        double t79 = -Ityy;
        double t83 = -Itzz;
        double t86 = Itzz*8.0;
        double t89 = Iwxx*8.0;
        double t90 = Iwyy*8.0;
        double t91 = -Iwzz;
        double t94 = Iwzz*8.0;
        double t97 = -Tneta;
        double t50 = cos(t41);
        double t51 = cos(t42);
        double t52 = cos(t43);
        double t53 = cos(t44);
        double t54 = t2*t2;
        double t55 = t3*t3;
        double t56 = t4*t4;
        double t57 = t4*t4*t4;
        double t58 = t5*t5;
        double t59 = t11*2.0;
        double t60 = t11*4.0;
        double t61 = t12*2.0;
        double t62 = sin(t41);
        double t63 = sin(t42);
        double t64 = sin(t43);
        double t65 = sin(t44);
        double t66 = t6*t6;
        double t67 = t7*t7;
        double t68 = t8*t8;
        double t69 = t9*t9;
        double t70 = mw+t10;
        double t71 = -t14;
        double t72 = -t15;
        double t76 = -t17;
        double t77 = -t18;
        double t80 = -t22;
        double t81 = -t23;
        double t82 = -t24;
        double t84 = -t25;
        double t85 = -t26;
        double t87 = -t28;
        double t88 = -t29;
        double t92 = -t31;
        double t93 = -t32;
        double t95 = rw*t11;
        double t96 = t6*dph;
        double t98 = rf*t10;
        double t99 = mt*t13;
        double t100 = t5*t6;
        double t101 = t40*2.0;
        double t102 = t6*t9;
        double t103 = t11*8.0;
        double t104 = t12*8.0;
        double t105 = rf*t11;
        double t106 = mt*t38;
        double t107 = mf*t40;
        double t108 = rt*t12;
        double t109 = mt*t40;
        double t110 = mw*t40;
        double t111 = t13*t13;
        double t112 = -t78;
        double t113 = -t86;
        double t114 = -t94;
        double t119 = t13*t34;
        double t120 = t13*t35;
        double t121 = -t48;
        double t123 = Itxx+t79;
        double t124 = Iwxx+t91;
        double t142 = t2*t5*t8;
        double t145 = t2*t8*t9;
        double t152 = t34+t36;
        double t153 = t35+t37;
        double t154 = Ifxx*t2*t4*t8;

        double t157 = Iwxx*t2*t3*t7;

        double t180 = t2*t4*t8*t75;
        double t181 = t2*t4*t8*t83;
        double t182 = t2*t3*t7*t91;
        double t188 = t10*2.0+t36;
        double t115 = rf*t61;
        double t116 = t96*2.0;
        double t117 = t96*4.0;
        double t118 = t98*2.0;
        double t122 = t46*t96;
        double t125 = rw*t70;
        double t126 = rw*t99;
        double t127 = t96+dps;
        double t128 = t96+dom;
        double t129 = t50*3.0;
        double t130 = rf*t59;
        double t131 = t105*3.0;
        double t132 = rf*t60;
        double t133 = -t96;
        double t134 = Itxx*t58;
        double t135 = Ityy*t58;
        double t136 = Iwxx*t55;
        double t137 = Ifzz*t68;
        double t138 = Itxx*t69;
        double t139 = Ityy*t69;
        double t140 = Itzz*t68;
        double t141 = Iwzz*t67;
        double t143 = t98*8.0;
        double t144 = t99*8.0;
        double t146 = t13*t99;
        double t149 = rf*t103;
        double t150 = t16*t66;
        double t151 = t30*t66;
        double t159 = t34*t111;
        double t161 = t35*t111;
        double t162 = t47*t50;
        double t163 = t47*t54;
        double t165 = t19+t80;
        double t166 = t20+t81;
        double t167 = t21+t82;
        double t168 = t123*t123;
        double t169 = t28+t92;
        double t170 = t29+t93;
        double t171 = t40*t70;
        double t172 = t45+t96;
        double t173 = t46+t96;
        double t174 = -t142;
        double t175 = t38+t101;
        double t177 = t12+t98;
        double t178 = t11+t99;
        double t187 = t89+t114;
        double t189 = t52*t54*2.0;
        double t191 = t18*t54*t56;
        double t192 = t26*t54*t56;
        double t193 = t32*t54*t55;
        double t194 = t40*t152;
        double t195 = t40*t153;
        double t196 = t15*t54*t68;
        double t197 = t29*t54*t67;
        double t200 = t59+t119;
        double t201 = t60+t120;
        double t202 = t53*t123;
        double t203 = t51*t124;
        double t204 = t63*t124;
        double t208 = t8*t53*t62*4.0;
        double t219 = t100+t145;
        double t224 = t6*t65*t123;
        double t230 = t5*t9*t49*t123;
        double t231 = t3*t7*t48*t124;
        double t234 = t40*t188;

        double t254 = t2*t4*t65*t123*2.0;
        double t255 = t62*t65*t123*dph;
        double t258 = t4*t62*t65*t123;
        double t259 = t2*t64*t65*t123;

        double t291 = t5*t9*t56*t121*t123;
        double t295 = t2*t5*t9*t56*t83*t123;
        double t147 = t117+dps;
        double t148 = -t129;
        double t160 = t146*3.0;

        double t176 = t13*t144;
        double t179 = t171*4.0;
        double t183 = t96*t127;
        double t184 = t46+t116+dps;
        double t185 = -t162;
        double t186 = Iwyy+t171;
        double t190 = t178*t178;
        double t198 = mf*t175;
        double t199 = t61+t118;
        double t205 = t194*4.0;
        double t206 = t33*t175;
        double t209 = t4*t178;
        double t210 = t202*2.0;
        double t211 = t202*4.0;
        double t212 = t203*4.0;
        double t213 = t53*t166;
        double t214 = rw*t4*t177;
        double t215 = t63*t169;
        double t217 = t105+t146;
        double t218 = Itzz+t202;
        double t220 = t104+t143;
        double t221 = t103+t144;

        double t225 = -t202;
        double t229 = t51*t169*2.0;
        double t233 = t224*4.0;

        double t236 = g*t2*t8*t178;
        double t240 = rw*t4*t200;
        double t241 = rw*t8*t200;
        double t242 = rw*t8*t201;
        double t243 = t130+t159;
        double t244 = rw*t201*dth*dom;
        double t246 = t102+t174;
        double t248 = t4*t224;
        double t250 = t219*t219;
        double t251 = rw*t8*t49*t177;
        double t253 = t54*t204*dps;
        double t257 = t2*t4*t65*t165;
        double t260 = t65*t96*t165*dth;
        double t262 = t9*t123*t174;

        double t264 = t63*t187*dth*dps;
        double t265 = t65*t123*t128;
        double t269 = t56*t65*t165*dth;
        double t270 = t259*2.0;
        double t275 = t3*t7*t124*t163;
        double t276 = t259*dph;
        double t280 = Ityy*t4*t9*t219;
        double t282 = rw*t8*t133*t177*dom;
        double t297 = t8*t50*t65*t167*dph;
        double t306 = t56*t58*t69*t168;
        double t316 = t2*t51*t170*t172*dph;
        double t327 = Itzz*rw*t2*t5*t9*t57*t123*t178;
        double t335 = Ifxx+t105+t106+t108+t115+t134+t139;
        double t336 = Ifyy+t105+t106+t108+t115+t135+t138;
        double t207 = t198*4.0;
        double t216 = rw*t209;
        double t223 = t30+t179;
        double t226 = -t210;

        double t228 = t213*2.0;
        double t232 = t48+t183;
        double t238 = -t213;
        double t245 = Itzz+t225;

        double t252 = t4*t217;
        double t256 = t25+t210;
        double t261 = -t236;

        double t268 = t246*t246;
        double t271 = -t251;
        double t272 = rw*t6*t8*t199*2.0;
        double t273 = fcoeff+t253;
        double t277 = t24*t250;
        double t278 = rw*t2*t220*dth;

        double t284 = t40*t56*t190;
        double t285 = t66*t242;
        double t286 = t2*t8*t218*4.0;
        double t287 = t125+t209;
        double t288 = t148+t189+1.0;
        double t289 = t6*t45*t241;
        double t290 = -t269;
        double t292 = -t275;
        double t296 = t49+t122+t185;
        double t305 = rw*t4*t66*t221;
        double t307 = Itxx*t4*t5*t246;
        double t308 = -t297;

        double t314 = t52*t54*t243;
        double t317 = rw*t2*t4*t147*t201;

        double t328 = rw*t2*t4*t184*t221*dph;
        double t329 = t146+t194+t198;
        double t331 = t186*t295;
        double t337 = t186+t203+t214;
        double t349 = t47*t123*t219*t246;
        double t353 = t56*t335;
        double t354 = Itxx+Ityy+t14+t76+t84+t202+t243;
        double t360 = t75+t83+t335;
        double t365 = t15+t19+t22+t77+t85+t132+t161+t210;
        double t369 = t186*t336;
        double t376 = t21+t24+t73+t112+t113+t149+t176+t211;
        double t239 = -t228;
        double t247 = t223*dps;
        double t267 = t25+t226;
        double t283 = t4*t245*dth;
        double t293 = t21*t268;
        double t294 = t2*t8*t245;
        double t298 = t4*t6*t256;
        double t301 = t273*dph*4.0;

        double t303 = rw*t287;
        double t313 = -t307;
        double t315 = t65*t288;
        double t318 = g*t6*t287*8.0;
        double t320 = t4*t128*t256;
        double t321 = -t314;

        double t323 = rw*t200*t232;
        double t324 = t95+t126+t252;
        double t326 = -dps*(t212-t223);
        double t332 = t65*t123*t296;

        double t339 = t254+t272;
        double t340 = t50*t329;
        double t341 = t255+t278;
        double t343 = t16+t19+t22+t132+t161+t226;

        double t355 = t2*t337*dph*dth;
        double t364 = t64*t354;
        double t371 = t64*t360;
        double t372 = t216+t336;
        double t381 = (t2*t52*t365*dph*dth)/4.0;
        double t382 = t64*t376*dth*dom;

        double t385 = t2*t52*t173*t365*dph;
        double t387 = t2*t52*t365*(t96-dom);
        double t403 = Itxx+Ityy+t16+t30+t71+t76+t84+t87+t92+t159+t195+t206+t238;

        double t415 = t107+t109+t110+t136+t137+t140+t141+t240+t353;
        double t299 = t294*dph;
        double t300 = -t283;
        double t304 = t4*t267*dom;
        double t309 = t4*t6*t267*2.0;
        double t310 = Iwyy+t303;
        double t333 = t2*t8*t324;
        double t334 = t224+t294;
        double t344 = -t340;

        double t350 = t343*dom;
        double t351 = t259+t298;
        double t352 = t49*t339;
        double t361 = (t4*t341*dph)/4.0;
        double t367 = (t2*t343*dph*dth)/4.0;
        double t368 = t260+t323;
        double t370 = -t123*dph*(t208-t315);
        double t373 = t244+t332;
        double t374 = t54*t364;

        double t377 = t6*t372;

        double t390 = (t364*(t48-t163))/4.0;
        double t392 = t276+t290+t320;

        double t401 = t215+t242+t364;
        double t402 = t204+t241+t371;
        double t408 = t96*t403;
        double t419 = t19+t22+t72+t74+t77+t85+t88+t90+t93+t161+t205+t207+t229+t239;
        double t423 = Itzz*t2*t4*t415;
        double t431 = t231+t271+t282+t292+t355+Tneta;
        double t433 = t284*t415;

        double t440 = t336*t415;

        double t319 = t310*t310;
        double t342 = t334*dph*dom;
        double t346 = t4*t9*t100*t123*t310;
        double t347 = t6*t216*t310;
        double t357 = t351*dph;
        double t358 = -t352;

        double t378 = -t374;
        double t379 = (t8*t368)/2.0;
        double t380 = t8*t373*4.0;
        double t384 = t265+t299+t300;
        double t386 = t6*t310*t336;
        double t395 = t392*dne*4.0;
        double t397 = t262+t377;
        double t404 = t6*t401;
        double t405 = t2*t402;
        double t424 = t96*t419;
        double t437 = t154+t157+t180+t181+t182+t280+t313+t333;

        double t448 = t247+t350+t408;
        double t449 = dne*(t370+dom*(t233-t286)+dth*(t270-t309));

        double t481 = t131+t150+t151+t160+t191+t192+t193+t196+t197+t234+t277+t293+t305+t321+t344;

        double t330 = t66*t319;
        double t356 = -t347;

        double t363 = t8*t83*t347;
        double t388 = t384*dne;

        double t394 = t8*t83*t386;
        double t398 = -t395;
        double t399 = dth*(t304-t357)*(-1.0/2.0);
        double t400 = Itzz*t8*t397;
        double t407 = t186*t397;
        double t409 = t4*t5*t9*t123*t397;
        double t410 = t216*t397;
        double t412 = t258+t285+t378;
        double t416 = t257+t404;
        double t417 = t6*t310*t397;

        double t426 = Itzz*t2*t4*(t248-t405)*(-1.0/2.0);
        double t429 = t4*t5*t9*t123*(t248-t405)*(-1.0/2.0);
        double t430 = t216*(t248-t405)*(-1.0/2.0);
        double t439 = Itzz*t8*t437;

        double t442 = t347*t415;
        double t444 = t336*(t248-t405)*(-1.0/2.0);
        double t446 = t4*t5*t9*t123*t437;
        double t447 = t216*t437;
        double t450 = t2*t448*dph*2.0;
        double t457 = t336*t437;
        double t458 = t397*t415;

        double t463 = t326+t350+t424;
        double t475 = t397*t437;

        double t482 = (t437*(t248-t405))/2.0;
        double t485 = (Itzz*t8*t481)/4.0;
        double t487 = (t186*t481)/4.0;
        double t489 = (t4*t5*t9*t123*t481)/4.0;
        double t490 = (t216*t481)/4.0;

        double t512 = (t415*t481)/4.0;
        double t411 = -t410;
        double t413 = t412*dph*2.0;
        double t414 = t295+t400;
        double t421 = t48*t416;

        double t438 = -Itzz*t8*(t347-t407);
        double t445 = t356*t415;
        double t454 = t346+t430;
        double t464 = t2*t463;
        double t466 = t346+t447;
        double t469 = t230+t291+t342+t349+t399+Tomega;

        double t478 = t409+t444;
        double t479 = ((t248-t405)*(t347-t407))/2.0;
        double t484 = t423+t439;
        double t488 = -t487;
        double t491 = -t490;
        double t492 = t409+t457;
        double t498 = t429+t458;
        double t508 = t446+t458;
        double t513 = t426+t485;
        double t516 = t97+t261+t361+t367+t379+t381+t388+t390;

        double t525 = t4*t5*t9*t123*(t489-(t397*(t248-t405))/2.0);
        double t527 = t475+t489;
        double t537 = t264+t316+t318+t328+t380+t382+t385+t398+t450;
        double t539 = t482+t512;

        double t422 = t8*t83*t414;
        double t425 = t186*t414;
        double t427 = t216*t414;
        double t428 = t289+t413;
        double t451 = t386+t411;
        double t456 = t331+t438;
        double t459 = t216*t454;
        double t470 = t4*t5*t9*t123*t466;
        double t471 = t216*t466;
        double t476 = t6*t310*t466;
        double t483 = t186*t478;
        double t486 = Itzz*t2*t4*t484;
        double t495 = t186*t492;
        double t496 = t6*t310*t492;
        double t497 = t330+t488;
        double t499 = t308+t317+t387+t464;
        double t504 = t186*t498;
        double t506 = t216*t498;
        double t509 = t6*t310*t498;
        double t511 = t186*t508;
        double t514 = t417+t491;
        double t515 = Itzz*t8*t513;
        double t520 = t397*t508;
        double t529 = t186*t527;
        double t531 = t216*t527;

        double t541 = t83*t539;
        double t542 = t216*t539;
        double t544 = t336*t539;
        double t432 = t428*dom;
        double t452 = t363+t425;
        double t453 = Itzz*t8*t451;

        double t465 = t394+t427;

        double t493 = t415*t451;

        double t500 = t4*t5*t9*t123*t497;
        double t502 = t499*dth;

        double t510 = t336*t497;
        double t517 = t4*t5*t9*t123*t514;
        double t518 = t216*t514;

        double t533 = -t531;
        double t534 = t445+t511;
        double t543 = -t542;
        double t545 = -t544;

        double t559 = t486+t515+t541;
        double t434 = -t432;
        double t467 = t327+t453;
        double t503 = -t500;

        double t549 = t496+t533;

        double t556 = t509+t543;
        double t560 = t216*t559;

        double t567 = Itzz*t8*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))*8.0;
        double t568 = t520+t525+t545;
        double t538 = t479+t503;
        double t548 = t301+t358+t421+t434+t449+t502;

        double t557 = t216*t556;

        double t569 = t186*t568;

        double t576 = -1.0/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double t577 = -1.0/(t26*(t557-t569+t6*t310*(t470-t493))+t8*t85*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+t2*t4*t26*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double t579 = -1.0/(-t567+t86*(t557-t569+t6*t310*(t470-t493))+t2*t4*t86*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));

        double et1 = (t548*(Itzz*(t433+t186*(t306-t440))-t27*t68*(t284-t369)))/(t26*(t557-t569+t6*t310*(t470-t493))+t8*t85*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+t2*t4*t26*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t431*t576*(Itzz*(t470-t493)-Itzz*t8*t465);
        double et2 = (t469*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t537*t579*(Itzz*(t471-t495)+t2*t4*t8*t27*(t284-t369));
        double et3 = (t516*(Itzz*t534+t8*t83*t452))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et4 = (t537*(Itzz*(-t510+t518+t397*(t347-t407))+t27*t54*t56*(t284-t369)))/(-t567+t86*(t557-t569+t6*t310*(t470-t493))+t2*t4*t86*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+(t516*(Itzz*t538+Itzz*t2*t4*t452))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et5 = (t431*(Itzz*(t517-(t451*(t248-t405))/2.0)+t2*t4*t83*t465))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t469*t576*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483));
        double et6 = (t548*(Itzz*(t459-t483)+t2*t4*t8*t27*(t284-t369)))/(t26*(t557-t569+t6*t310*(t470-t493))+t8*t85*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+t2*t4*t26*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et7 = (t431*(Itzz*t568+t8*t83*(t397*t414-t336*t513)+Itzz*t2*t4*(t336*t484+t4*t5*t9*t123*t414)))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t516*t576*(t560+t6*t310*(t422+Itzz*t508))+t537*t579*(Itzz*t549+Itzz*t2*t4*t467)+t548*t577*(Itzz*(t506+t6*t310*(t306-t440))-t8*t83*t467);
        double et8 = (t469*(Itzz*t8*t549-Itzz*t2*t4*(t506+t6*t310*(t306-t440))))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et9 = (t516*(t186*t559+t330*(Itzz*t107+Itzz*t109+Itzz*t110+Itzz*t136+Itzz*t137+Itzz*t141+Ifxx*Itzz*t56+Itzz*t56*t105+Itzz*t56*t106+Itzz*t56*t108+Itzz*t56*t134+Itzz*t56*t139+t4*t25*t95+rf*t12*t25*t56+rw*t4*t12*t25+mt*rf*rw*t4*t25)))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))));
        double et10 = (t537*(Itzz*(t476-t529)-Itzz*t2*t4*t456))/(-t567+t86*(t557-t569+t6*t310*(t470-t493))+t2*t4*t86*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t431*t576*(t560+t6*t310*(t422+Itzz*t498))+t469*t576*(Itzz*t8*(t476-t529)+Itzz*t2*t4*(t442-t504))+t548*t577*(Itzz*(t442-t504)-t8*t83*t456);
        double et11 = t516*t576*(Itzz*t8*t538+Itzz*t2*t4*t534)+(t469*(t557-t569+t6*t310*(t470-t493)))/(Itzz*(t557-t569+t6*t310*(t470-t493))+t8*t83*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t459-t483))+Itzz*t2*t4*(Itzz*t8*(t471-t495)+Itzz*t2*t4*(t433+t186*(t306-t440))))+t537*t579*(Itzz*t8*(-t510+t518+t397*(t347-t407))+t2*t4*t83*(t471-t495))+t548*t577*(Itzz*t8*(t459-t483)+Itzz*t2*t4*(t433+t186*(t306-t440)));
        double et12 = t431*t576*(Itzz*t8*(t517-(t451*(t248-t405))/2.0)-Itzz*t2*t4*(t470-t493));

        d_dNE[index] = et11+et12;
        
        index = index + num_threads;
        index1 = (grid_size[0] > 1) ? index : 0;
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
        } else if (curr_free_dim == 1) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_THn, p_TH, p_k1_TH, p_k2_TH, p_k3_TH, p_k4_TH, dt, p_limits[1], p_limits[9], p_grid_size);
            mxGPUDestroyGPUArray(k1_TH);
            mxGPUDestroyGPUArray(k2_TH);
            mxGPUDestroyGPUArray(k3_TH);
            mxGPUDestroyGPUArray(k4_TH);
        } else if (curr_free_dim == 2) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_OMn, p_OM, p_k1_OM, p_k2_OM, p_k3_OM, p_k4_OM, dt, p_limits[2], p_limits[10], p_grid_size);
            mxGPUDestroyGPUArray(k1_OM);
            mxGPUDestroyGPUArray(k2_OM);
            mxGPUDestroyGPUArray(k3_OM);
            mxGPUDestroyGPUArray(k4_OM);
        } else if (curr_free_dim == 3) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_dPHn, p_dPH, p_k1_dPH, p_k2_dPH, p_k3_dPH, p_k4_dPH, dt, p_limits[3], p_limits[11], p_grid_size);
            mxGPUDestroyGPUArray(k1_dPH);
            mxGPUDestroyGPUArray(k2_dPH);
            mxGPUDestroyGPUArray(k3_dPH);
            mxGPUDestroyGPUArray(k4_dPH);
        } else if (curr_free_dim == 4) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_dTHn, p_dTH, p_k1_dTH, p_k2_dTH, p_k3_dTH, p_k4_dTH, dt, p_limits[4], p_limits[12], p_grid_size);
            mxGPUDestroyGPUArray(k1_dTH);
            mxGPUDestroyGPUArray(k2_dTH);
            mxGPUDestroyGPUArray(k3_dTH);
            mxGPUDestroyGPUArray(k4_dTH);
        } else if (curr_free_dim == 5) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_dPSn, p_dPS, p_k1_dPS, p_k2_dPS, p_k3_dPS, p_k4_dPS, dt, p_limits[5], p_limits[13], p_grid_size);
            mxGPUDestroyGPUArray(k1_dPS);
            mxGPUDestroyGPUArray(k2_dPS);
            mxGPUDestroyGPUArray(k3_dPS);
            mxGPUDestroyGPUArray(k4_dPS);
        } else if (curr_free_dim == 6) {
            rk4_wbounds << <blocksPerGrid, threadsPerBlock >> > (p_dOMn, p_dOM, p_k1_dOM, p_k2_dOM, p_k3_dOM, p_k4_dOM, dt, p_limits[6], p_limits[14], p_grid_size);
            mxGPUDestroyGPUArray(k1_dOM);
            mxGPUDestroyGPUArray(k2_dOM);
            mxGPUDestroyGPUArray(k3_dOM);
            mxGPUDestroyGPUArray(k4_dOM);
        } else if (curr_free_dim == 7) {
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
