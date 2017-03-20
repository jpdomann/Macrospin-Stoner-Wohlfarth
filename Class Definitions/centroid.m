function [cg,area]=centroid(p)

    % make a copy of the pts offset by 1
    q=p([2:end 1],:);

    % compute partial terms
    da = p(:,1).*q(:,2) - q(:,1).*p(:,2);
    dx = (q(:,1) + p(:,1)) .* da;
    dy = (q(:,2) + p(:,2)) .* da;

    % sum
    a = sum(da);
    x = sum(dx);
    y = sum(dy);

    % use those to compute signed area & centroid
    area = a / 2;
    cg = [x/(3*a),y/(3*a)];
end