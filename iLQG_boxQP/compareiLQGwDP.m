clear;
clc;
close all;

%% 
iLQG_dir = '';
iLQG_filename = 'iLQGDecomposed_test.mat';
% DP_dir = 'data/';
% DP_filename = 'mc=5,mp=1_densegrid.mat';
DP_cont_dir = 'data/';
DP_cont_filename = 'mc=5,mp=1_continuousa_densegrid.mat';

load(strcat(iLQG_dir, iLQG_filename), 'XJoint', 'XCartFS', 'XPoleTS', 'XCartTS', 'XPoleFS', ...
                                      'XCartFPoleTDec', 'XCartTPoleFDec', ...
                                      'XCartFClose', 'XPoleTClose', 'XCartTClose', 'XPoleFClose', ...
                                      'UCartFClose', 'UPoleTClose', 'UCartTClose', 'UPoleFClose', ...
                                      'KCartFClose', 'KPoleTClose', 'KCartTClose', 'KPoleFClose', ...
                                      'UJoint', 'UCartFS', 'UPoleTS', 'UCartTS', 'UPoleFS', ...
                                      'UCartFPoleTDec', 'UCartTPoleFDec', 'sys', 'x_starts');
load(strcat(DP_cont_dir, DP_cont_filename), 'V_joint', 'policy_joint', ...
                                  'V_t_pole_second', 'f_cart_first_', 't_pole_second', ...
                                  'V_f_cart_second', 't_pole_first_', 'f_cart_second', ...
                                  'V_f_pole_second', 't_cart_first_', 'f_pole_second', ...
                                  'V_t_cart_second', 'f_pole_first_', 't_cart_second', ...
                                  'grid_x', 'grid_x_dot', 'grid_xP', 'grid_xP_dot', ...
                                  'x_limits', 'x_dot_limits', 'xP_limits', 'xP_dot_limits');

%% Estimate DDP values
J_DDP = nan(size(x_starts,2), 7);

for jj=1:1:size(x_starts, 2)
    
    J_DDP(jj, 1) = J_CartPole_WrapAround(sys, XJoint(:,:, jj), UJoint(:,:, jj));
    
    sys.U_DIMS_FREE = [2];
    cost         = l_CartPoleSecondWrapAround(sys, XPoleTS(:,:, jj), UPoleTS(:,:, jj), ...
                                              UCartFClose(:,:, jj), KCartFClose(:,:,:, jj), XCartFClose(:,:, jj));
    discount     = sys.gamma_.^linspace(0, length(cost)-1, length(cost));
    J_DDP(jj, 2) = sum(cost.*discount*sys.dt, 2);
    
    sys.U_DIMS_FREE = [1];
    cost         = l_CartPoleSecondWrapAround(sys, XCartFS(:,:, jj), UCartFS(:,:, jj), ...
                                              UPoleTClose(:,:, jj), KPoleTClose(:,:,:, jj), XPoleTClose(:,:, jj));
    discount     = sys.gamma_.^linspace(0, length(cost)-1, length(cost));
    J_DDP(jj, 3) = sum(cost.*discount*sys.dt, 2);
    
    sys.U_DIMS_FREE = [1];
    cost         = l_CartPoleSecondWrapAround(sys, XPoleFS(:,:, jj), UPoleFS(:,:, jj), ...
                                              UCartTClose(:,:, jj), KCartTClose(:,:,:, jj), XCartTClose(:,:, jj));
    discount     = sys.gamma_.^linspace(0, length(cost)-1, length(cost));
    J_DDP(jj, 4) = sum(cost.*discount*sys.dt, 2);
    
    sys.U_DIMS_FREE = [2];
    cost         = l_CartPoleSecondWrapAround(sys, XCartTS(:,:, jj), UCartTS(:,:, jj), ...
                                              UPoleFClose(:,:, jj), KPoleFClose(:,:,:, jj), XPoleFClose(:,:, jj));
    discount     = sys.gamma_.^linspace(0, length(cost)-1, length(cost));
    J_DDP(jj, 5) = sum(cost.*discount*sys.dt, 2);
    
    J_DDP(jj, 6) = J_CartPole_WrapAround(sys, XCartFPoleTDec(:,:, jj), UCartFPoleTDec(:,:, jj));
    J_DDP(jj, 7) = J_CartPole_WrapAround(sys, XCartTPoleFDec(:,:, jj), UCartTPoleFDec(:,:, jj));
end


%% DP
sys.U_DIMS_FREE = [1;2];
sys.X_DIMS_FREE = [1;2;3;4];
sys.dynamics_discrete = @(x, u) cartpole_dyn_first_cst(sys, x, u, sys.full_DDP);
limits = [x_limits; x_dot_limits; xP_limits; xP_dot_limits];

[trajXJoint, trajUJoint, JJoint] = rolloutDPTraj(policy_joint, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot);

policyCartFPoleT(:,:,:,:,1) = f_cart_first_;
policyCartFPoleT(:,:,:,:,2) = t_pole_second;
[trajXCFPT, trajUCFPT, JCFPT] = rolloutDPTraj(policyCartFPoleT, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot);

policyPoleTCartF(:,:,:,:,2) = t_pole_first_;
policyPoleTCartF(:,:,:,:,1) = f_cart_second;
[trajXPTCF, trajUPTCF, JPTCF] = rolloutDPTraj(policyPoleTCartF, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot);

policyCartTPoleF(:,:,:,:,2) = t_cart_first_;
policyCartTPoleF(:,:,:,:,1) = f_pole_second;
[trajXCTPF, trajUCTPF, JCTPF] = rolloutDPTraj(policyCartTPoleF, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot);

policyPoleFCartT(:,:,:,:,1) = f_pole_first_;
policyPoleFCartT(:,:,:,:,2) = t_cart_second;
[trajXPFCT, trajUPFCT, JPFCT] = rolloutDPTraj(policyPoleFCartT, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot);

policyCartFPoleTDec(:,:,:,:,1) = f_cart_first_;
policyCartFPoleTDec(:,:,:,:,2) = t_pole_first_;
[trajXCFPTD, trajUCFPTD, JCFPTD] = rolloutDPTraj(policyCartFPoleTDec, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot);

policyCartTPoleFDec(:,:,:,:,1) = f_pole_first_;
policyCartTPoleFDec(:,:,:,:,2) = t_cart_first_;
[trajXCTPFD, trajUCTPFD, JCTPFD] = rolloutDPTraj(policyCartTPoleFDec, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot);

%% Compare trajectories

for jj=1:1:size(x_starts, 2)

    figure;
    subplot(2,1,1);
    hold on;
    % plot DDP trajectories
    plot(XJoint(1,:,jj),XJoint(2,:,jj),'Color',[0.4660 0.6740 0.1880]);
    plot(XPoleTS(1,:,jj),XPoleTS(2,:,jj),'Color',[0 0.4470 0.7410]);
    plot(XCartFS(1,:,jj),XCartFS(2,:,jj),'Color',[0.8500 0.3250 0.0980]);
    plot(XPoleFS(1,:,jj),XPoleFS(2,:,jj),'Color',[0.9290 0.6940 0.1250]);
    plot(XCartTS(1,:,jj),XCartTS(2,:,jj),'Color',[0.4940 0.1840 0.5560]);
    plot(XCartFPoleTDec(1,:,jj), XCartFPoleTDec(2,:,jj), 'Color', [0.3010 0.7450 0.9330]);
    plot(XCartTPoleFDec(1,:,jj), XCartTPoleFDec(2,:,jj), 'Color', [0.6350 0.0780 0.1840]);
    % plot DP trajectories
    plot(trajXJoint(1,:,jj), trajXJoint(2,:,jj),'Color',[0.4660 0.6740 0.1880], 'LineStyle', '-.');
    plot(trajXCFPT(1,:,jj), trajXCFPT(2,:,jj),'Color',[0 0.4470 0.7410], 'LineStyle', '-.');
    plot(trajXPTCF(1,:,jj), trajXPTCF(2,:,jj),'Color',[0.8500 0.3250 0.0980], 'LineStyle', '-.');
    plot(trajXCTPF(1,:,jj), trajXCTPF(2,:,jj),'Color',[0.9290 0.6940 0.1250], 'LineStyle', '-.');
    plot(trajXPFCT(1,:,jj), trajXPFCT(2,:,jj),'Color',[0.4940 0.1840 0.5560], 'LineStyle', '-.');
    plot(trajXCFPTD(1,:,jj), trajXCFPTD(2,:,jj), 'Color', [0.3010 0.7450 0.9330], 'LineStyle', '-.');
    plot(trajXCTPFD(1,:,jj), trajXCTPFD(2,:,jj), 'Color', [0.6350 0.0780 0.1840], 'LineStyle', '-.');
    
    legend(["Joint (DDP)", "F - Cart, T - Both (DDP)", "T - Pole, F - Both (DDP)", "T - Cart, F - Both (DDP)", "F - Pole, T - Both (DDP)", ...
            "F - Cart, T - Pole (DDP)", "T - Cart, F - Pole (DDP)", ...
            "Joint (DP)", "F - Cart, T - Both (DP)", "T - Pole, F - Both (DP)", "T - Cart, F - Both (DP)", "F - Pole, T - Both (DP)", ...
            "F - Cart, T - Pole (DP)", "T - Cart, F - Pole (DP)"]);
    scatter(XJoint(1,1,jj),XJoint(2,1,jj), 20, [1,0,0]);
    xlabel('x');
    ylabel('x-dot');
    hold off;
    
    subplot(2,1,2);
    hold on;
    % plot DDP trajectories
    plot(XJoint(3,:,jj), XJoint(4,:,jj),'Color',[0.4660 0.6740 0.1880]);
    plot(XPoleTS(3,:,jj), XPoleTS(4,:,jj),'Color',[0 0.4470 0.7410]);
    plot(XCartFS(3,:,jj), XCartFS(4,:,jj),'Color',[0.8500 0.3250 0.0980]);
    plot(XPoleFS(3,:,jj), XPoleFS(4,:,jj),'Color',[0.9290 0.6940 0.1250]);
    plot(XCartTS(3,:,jj), XCartTS(4,:,jj),'Color',[0.4940 0.1840 0.5560]);
    plot(XCartFPoleTDec(3,:,jj), XCartFPoleTDec(4,:,jj), 'Color', [0.3010 0.7450 0.9330]);
    plot(XCartTPoleFDec(3,:,jj), XCartTPoleFDec(4,:,jj), 'Color', [0.6350 0.0780 0.1840]);
    % plot DP trajectories
    plot(trajXJoint(3,:,jj), trajXJoint(4,:,jj),'Color',[0.4660 0.6740 0.1880], 'LineStyle', '-.');
    plot(trajXCFPT(3,:,jj), trajXCFPT(4,:,jj),'Color',[0 0.4470 0.7410], 'LineStyle', '-.');
    plot(trajXPTCF(3,:,jj), trajXPTCF(4,:,jj),'Color',[0.8500 0.3250 0.0980], 'LineStyle', '-.');
    plot(trajXCTPF(3,:,jj), trajXCTPF(4,:,jj),'Color',[0.9290 0.6940 0.1250], 'LineStyle', '-.');
    plot(trajXPFCT(3,:,jj), trajXPFCT(4,:,jj),'Color',[0.4940 0.1840 0.5560], 'LineStyle', '-.');
    plot(trajXCFPTD(3,:,jj), trajXCFPTD(4,:,jj), 'Color', [0.3010 0.7450 0.9330], 'LineStyle', '-.');
    plot(trajXCTPFD(3,:,jj), trajXCTPFD(4,:,jj), 'Color', [0.6350 0.0780 0.1840], 'LineStyle', '-.');
    
    legend(["Joint (DDP)", "F - Cart, T - Both (DDP)", "T - Pole, F - Both (DDP)", "T - Cart, F - Both (DDP)", "F - Pole, T - Both (DDP)", ...
            "F - Cart, T - Pole (DDP)", "T - Cart, F - Pole (DDP)", ...
            "Joint (DP)", "F - Cart, T - Both (DP)", "T - Pole, F - Both (DP)", "T - Cart, F - Both (DP)", "F - Pole, T - Both (DP)", ...
            "F - Cart, T - Pole (DP)", "T - Cart, F - Pole (DP)"]);
    scatter(XJoint(3,1,jj),XJoint(4,1,jj),20,[1,0,0]);
    xlabel('theta'); xlims()
    ylabel('theta-dot');
    hold off;
end

function [trajX, trajU, J] = rolloutDPTraj(policy, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot)
    
    NUM_CTRL = round(sys.T / sys.dt);
    trajX = nan(4, NUM_CTRL+1, size(x_starts, 2));
    trajU = nan(2, NUM_CTRL, size(x_starts, 2));
    J = zeros(size(x_starts, 2), 1);
    
    for jj = 1:1:size(x_starts, 2)
        trajX(:, 1, jj) = x_starts(: ,jj);
        for ii=1:1:NUM_CTRL
            x_ = [min(limits(1, 2), max(limits(1, 1), trajX(1, ii, jj)));
                  min(limits(2, 2), max(limits(2, 1), trajX(2, ii, jj)));
                  trajX(3, ii, jj);
                  min(limits(4, 2), max(limits(4, 1), trajX(4, ii, jj)))];
            while((x_(3) > limits(3, 2)))
                x_(3) = x_(3) - limits(3, 2);
            end
            while ((x_(3) < limits(3, 1)))
                x_(3) = x_(3) + limits(3, 2);
            end
            trajU(1, ii, jj) = interpn(grid_x, grid_x_dot, grid_xP, grid_xP_dot, policy(:,:,:,:,1), ...
                                       x_(1),  x_(2),  x_(3),  x_(4));
            trajU(2, ii, jj) = interpn(grid_x, grid_x_dot, grid_xP, grid_xP_dot, policy(:,:,:,:,2), ...
                                       x_(1),  x_(2),  x_(3),  x_(4));
            [trajX(:, ii+1, jj),~] = sys.dynamics_discrete(trajX(:, ii, jj), trajU(:, ii, jj));
        end
        J(jj, 1) = J_CartPole_WrapAround(sys, trajX(:,:, jj), trajU(:,:, jj));
    end
end

