function [p_joints] = fcn_p_joints(q,p)

p_joints = zeros(2,6);

  p_joints(1,1)=q(1);
  p_joints(1,2)=q(1) + (p(2)*sin(q(3)))/2;
  p_joints(1,3)=q(1) + p(3)*sin(q(4) + q(3)) + (p(2)*sin(q(3)))/2;
  p_joints(1,4)=q(1) + p(3)*sin(q(4) + q(3)) + (p(2)*sin(q(3)))/2 + p(4)*sin(q(4) + q(5) + q(3));
  p_joints(1,5)=q(1) + p(3)*sin(q(6) + q(3)) + (p(2)*sin(q(3)))/2;
  p_joints(1,6)=q(1) + p(3)*sin(q(6) + q(3)) + (p(2)*sin(q(3)))/2 + p(4)*sin(q(6) + q(7) + q(3));
  p_joints(2,1)=q(2);
  p_joints(2,2)=q(2) - (p(2)*cos(q(3)))/2;
  p_joints(2,3)=q(2) - p(3)*cos(q(4) + q(3)) - (p(2)*cos(q(3)))/2;
  p_joints(2,4)=q(2) - p(3)*cos(q(4) + q(3)) - (p(2)*cos(q(3)))/2 - p(4)*cos(q(4) + q(5) + q(3));
  p_joints(2,5)=q(2) - p(3)*cos(q(6) + q(3)) - (p(2)*cos(q(3)))/2;
  p_joints(2,6)=q(2) - p(3)*cos(q(6) + q(3)) - (p(2)*cos(q(3)))/2 - p(4)*cos(q(6) + q(7) + q(3));

 