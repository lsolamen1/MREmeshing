function [in]=Hex27incidencelist(nodnum)
%Hex27incidencelist:outputs the row of the incidence list from a Hex27
%element built of of a 3x3x3 cube of node numbers, nodnum.

% Nodal coords copied directly from matts masters thesis.
% X(1,:) = [-1,-1,-1 ]; X(10,:) = [-1,-1, 1 ]; X(19,:) = [-1,-1, 0 ];
% X(2,:) = [ 1,-1,-1 ]; X(11,:) = [ 1,-1, 1 ]; X(20,:) = [ 1,-1, 0 ];
% X(3,:) = [ 1, 1,-1 ]; X(12,:) = [ 1, 1, 1 ]; X(21,:) = [ 1, 1, 0 ];
% X(4,:) = [-1, 1,-1 ]; X(13,:) = [-1, 1, 1 ]; X(22,:) = [-1, 1, 0 ];
% X(5,:) = [ 0,-1,-1 ]; X(14,:) = [ 0,-1, 1 ]; X(23,:) = [ 0,-1, 0 ];
% X(6,:) = [ 1, 0,-1 ]; X(15,:) = [ 1, 0, 1 ]; X(24,:) = [ 1, 0, 0 ];
% X(7,:) = [ 0, 1,-1 ]; X(16,:) = [ 0, 1, 1 ]; X(25,:) = [ 0, 1, 0 ];
% X(8,:) = [-1, 0,-1 ]; X(17,:) = [-1, 0, 1 ]; X(26,:) = [-1, 0, 0 ];
% X(9,:) = [ 0, 0,-1 ]; X(18,:) = [ 0, 0, 1 ]; X(27,:) = [ 0, 0, 0 ];
% Plot the nodal coords to check
% for ii=1:27
%     plot3(X(ii,1),X(ii,2),X(ii,3),'r.','markersize',18)
%     title(['node ' int2str(ii) ' added'])
%     grid on
%     hold on
%     pause
% end
% Added 2 to X in matlab to get indices for nodnum

in=zeros(1,27);

in(1)=nodnum(1,1,1);
in(2)=nodnum(3,1,1);
in(3)=nodnum(3,3,1);
in(4)=nodnum(1,3,1);
in(5)=nodnum(2,1,1);
in(6)=nodnum(3,2,1);
in(7)=nodnum(2,3,1);
in(8)=nodnum(1,2,1);
in(9)=nodnum(2,2,1);
in(10)=nodnum(1,1,3);
in(11)=nodnum(3,1,3);
in(12)=nodnum(3,3,3);
in(13)=nodnum(1,3,3);
in(14)=nodnum(2,1,3);
in(15)=nodnum(3,2,3);
in(16)=nodnum(2,3,3);
in(17)=nodnum(1,2,3);
in(18)=nodnum(2,2,3);
in(19)=nodnum(1,1,2);
in(20)=nodnum(3,1,2);
in(21)=nodnum(3,3,2);
in(22)=nodnum(1,3,2);
in(23)=nodnum(2,1,2);
in(24)=nodnum(3,2,2);
in(25)=nodnum(2,3,2);
in(26)=nodnum(1,2,2);
in(27)=nodnum(2,2,2);
end