function J = J_CartPole_WrapAround(sys, X, U )
%CALC_J 
NUM_CTRL = size(X, 2) - 1;
U_DIM = size(U, 1);

J = 0;
discount = 1;
for t = 1 : NUM_CTRL
    J = J + discount*l_CartPole_WrapAround(sys, X(:, t), U(:, t))*sys.dt;
    discount = discount*sys.gamma_;
end
J = J + discount*l_CartPole_WrapAround(sys, X(:, NUM_CTRL+1), zeros(U_DIM, 1))*sys.dt;
end

