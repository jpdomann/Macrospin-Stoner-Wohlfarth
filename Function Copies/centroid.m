function [cg,area]=centroid(x1,y1)

    % make a copy of the pts offset by 1
    x2=x1([2:end 1]);
    y2=y1([2:end 1]);

    % compute partial terms
    da = x1.*y2 - x2.*y1;
    dx = (x2 + x1) .* da;
    dy = (y2 + y1) .* da;

    % sum
    a = sum(da);
    x = sum(dx);
    y = sum(dy);

    % use those to compute signed area & centroid
    area = a / 2;
    cg = [x/(3*a),y/(3*a)];
          
end