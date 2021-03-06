% This is a first attempt to make a 3D model from the data given to us by
% Will Dickson

% The joint of the fly wing is located at [0 -1.3412 0.1442]

XTrans = [0 -1.3412 0.1442];
load flywing
load flybody

body = flybody;
wing = flywing;
wing.pts = wing.pts - repmat(XTrans,size(wing.pts,1),1);

%Approximate translation from center of fly body to the wing joint location
%JTrans = [0 +/-0.35 0.5]
RJTrans = [-0.18 -0.36 0.52];
LJTrans = [-0.18 0.36 0.52];

bodyscale = 1;
wingscale = 1;

%Body
body.pts = bodyscale.*body.pts;

%Left Wing
Lwing = wing;
Lwing.pts = wingscale.*Lwing.pts;
Lwing.pts = Lwing.pts + repmat(bodyscale.*LJTrans,size(Lwing.pts,1),1);

%Right Wing
%Rotate about z-axis by 180 degrees
Rz = @(theta) [cos(theta) -sin(theta) 0
                sin(theta) cos(theta) 0
                0 0 1];
Rwing = wing;
Rwing.pts = wingscale.*Rwing.pts;
Rwing.pts = (Rz(pi)*Rwing.pts')' + repmat(bodyscale.*RJTrans,size(Rwing.pts,1),1);

%figure; 
% tt = 3;
% trisurf(body.TRI,body.pts(:,1),body.pts(:,2),body.pts(:,3)+tt,'facecolor','g');
% hold on;
% trisurf(Rwing.TRI,Rwing.pts(:,1),Rwing.pts(:,2),Rwing.pts(:,3)+tt,'facecolor','b');
% trisurf(Lwing.TRI,Lwing.pts(:,1),Lwing.pts(:,2),Lwing.pts(:,3)+tt,'facecolor','b');
% axis equal

tt = 2;
hold on;
trisurf(body.TRI,body.pts(:,3),body.pts(:,2)+tt,-body.pts(:,1),'facecolor','g');
trisurf(Rwing.TRI,Rwing.pts(:,3),Rwing.pts(:,2)+tt,-Rwing.pts(:,1),'facecolor','b');
trisurf(Lwing.TRI,Lwing.pts(:,3),Lwing.pts(:,2)+tt,-Lwing.pts(:,1),'facecolor','b');
axis equal


% %=========================
% %Grab planar points of wing.
% idx = find(flywing.pts(:,1) > 0);
% xx = flywing.pts(idx,2:3);
% %reorder points
% idx = [1 2 4 5 7 11 13 15 18 20 22 25 27 29 31 33 35 36 38 40 41 39 37 34 ... 
%     32 30 28 26 24 23 21 19 17 16 14 12 10 9 8 6 3];
% xx = xx(idx,:);
%  
%     
% %Fit a Closed spline to the boundary
% dd = reshape(xx,[],1);
% 
% % -- the number of control points be about one sixth of the total
% % data points
% c = 4; %cubic spline
% numrep = c-1; %this is the number of control points that are repeated
% %at the beginning and end of the vector;
% npts = max(c+numrep,round(size(xx,1)/2));
% %Parameterize the points by chord length;
% Dist = sqrt(sum(diff(xx,1,1).^2,2));
% cDist = cumsum(Dist)';
% T = [0 cDist]./cDist(end) - 0.5;
% knotv = oknotP(npts,c,1,0);
% [N,D1,D2] = dbasisP(c,T,npts,knotv);
% 
% Aeq = zeros(2*numrep,2*npts);
% for i = 1:2
%     Aeq((i-1)*numrep+(1:numrep),(i-1)*npts+(1:numrep)) = eye(numrep);
%     Aeq((i-1)*numrep+(1:numrep),(i-1)*npts+(npts-(numrep-1):npts)) = -1*eye(numrep);
% end
% options = optimset('largescale','off');
% bb = lsqlin(blkdiag(N,N),dd,[],[],Aeq,zeros(2*numrep,1),[],[],[],options);
% 
% B = reshape(bb,[],2);
% 
% T = linspace(-0.5,0.5,size(xx,1));
% [N,D1,D2] = dbasisP(c,T,npts,knotv);
% bndytemp = N*B;
% flygenmod.wing = B;
% 
% %
% %
% %==========================================================================
% %Grab X-Z plane points of body.
% %
% %
% %
% xzpts = body.pts((body.pts(:,2) > -1e-4 & body.pts(:,2) < 1e-4),[1 3]);
% 
% %Reorder the points;
% startidx = 87;
% allpts = xzpts;
% count = 1;
% xx = zeros(size(xzpts));
% xx(count,:) = allpts(startidx,:);
% allpts(startidx,:) = [];
% while ~isempty(allpts)
%     cpt_idx = kdtreeidx(allpts,xx(count,:));
%     xx(count+1,:) = allpts(cpt_idx,:);
%     allpts(cpt_idx,:) = [];
%     count = count+1;
% end
% 
% figure; plot(xx(:,1),xx(:,2),'*')
% axis equal
% for k = 1:size(xx,1)
% hold on; text(xx(k,1),xx(k,2),num2str(k));
% end
% % Fix the high curvature at the dorsal neck region.
% xx = [xx(1:92,:)
%     xx([283 282],:)
%     xx(93:279,:)]; 
% 
% %Interpolate in the tail region where there are no points
% newx = interp1([0 1],[xx(end,1);xx(1,1)],0.1:0.1:0.9)';
% newy = interp1([0 1],[xx(end,2);xx(1,2)],0.1:0.1:0.9)';
% 
% xx = [xx
%     newx newy];
% 
% %Fit a Closed spline to the boundary
% dd = reshape(xx,[],1);
% 
% % -- the number of control points be about one sixth of the total
% % data points
% c = 4; %cubic spline
% numrep = c-1; %this is the number of control points that are repeated
% %at the beginning and end of the vector;
% npts = max(c+numrep,round(size(xx,1)/6));
% %Parameterize the points by chord length;
% Dist = sqrt(sum(diff(xx,1,1).^2,2));
% cDist = cumsum(Dist)';
% T = [0 cDist]./cDist(end) - 0.5;
% %T = linspace(-0.5,0.5,size(xx,1));
% knotv = oknotP(npts,c,1,0);
% [N,D1,D2] = dbasisP(c,T,npts,knotv);
% 
% Aeq = zeros(2*numrep,2*npts);
% for i = 1:2
%     Aeq((i-1)*numrep+(1:numrep),(i-1)*npts+(1:numrep)) = eye(numrep);
%     Aeq((i-1)*numrep+(1:numrep),(i-1)*npts+(npts-(numrep-1):npts)) = -1*eye(numrep);
% end
% options = optimset('largescale','off');
% bb = lsqlin(blkdiag(N,N),dd,[],[],Aeq,zeros(2*numrep,1),[],[],[],options);
% 
% B = reshape(bb,[],2);
% 
% T = linspace(-0.5,0.5,size(xx,1));
% [N,D1,D2] = dbasisP(c,T,npts,knotv);
% bndytemp = N*B;
% 
% flygenmod.bodybndy = B;


