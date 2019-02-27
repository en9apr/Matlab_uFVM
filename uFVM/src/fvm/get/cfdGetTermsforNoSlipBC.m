function [e_x,e_y,mag_Sb,d_perp]=cfdgettermsforNoslipBC(int_pts,centroid,iInterfaceArray)


% Distance between each liq/sol centroid point and corresponding piecewise linearized curve-line and finding the equation for the normal
% line


for i=1:length(iInterfaceArray)
%  hold on;
%  plot([int_pts(2*i-1,1)' int_pts(2*i,1)' centroid(iInterfaceArray(i),1) int_pts(2*i-1,1)'],...
%      [int_pts(2*i-1,2)' int_pts(2*i,2)' centroid(iInterfaceArray(i),2) int_pts(2*i-1,2)']);
 
 
 mag_Sb(i)=sqrt((int_pts(2*i-1,1)- int_pts(2*i,1))^2+(int_pts(2*i-1,2)- int_pts(2*i,2))^2);
 
 u(i)=((centroid(iInterfaceArray(i),1)-int_pts(2*i-1,1))*(int_pts(2*i,1)-int_pts(2*i-1,1)) ...
      +(centroid(iInterfaceArray(i),2)-int_pts(2*i-1,2))*(int_pts(2*i,2)-int_pts(2*i-1,2)))/((int_pts(2*i-1,1)-int_pts(2*i,1))^2 ...
      +(int_pts(2*i-1,2)-int_pts(2*i,2))^2);
 



Normal_pt(i,:)=[int_pts(2*i-1,1)+u(i)*(int_pts(2*i,1)-int_pts(2*i-1,1)),int_pts(2*i-1,2)+u(i)*(int_pts(2*i,2)-int_pts(2*i-1,2))];


d_perp(i)=sqrt((Normal_pt(i,1)-centroid(iInterfaceArray(i),1))^2+(Normal_pt(i,2)-centroid(iInterfaceArray(i),2))^2);



% Direction of normal line

dy(i)=Normal_pt(i,2)-centroid(iInterfaceArray(i),2);
dx(i)=Normal_pt(i,1)-centroid(iInterfaceArray(i),1);


% Direction of normal line unit vector

e_y(i)=dy(i)/d_perp(i);
e_x(i)=dx(i)/d_perp(i);







% hold on;
%  plot([Normal_pt(i,1),Normal_pt(i,1)+e_x(i)],[Normal_pt(i,2),Normal_pt(i,2)+e_y(i)]);

end




% 
% Normal_pt(:,1)
% 
% sqrt((e_x).^2+(e_y).^2)



end
