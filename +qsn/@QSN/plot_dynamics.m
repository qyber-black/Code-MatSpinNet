function plot_dynamics (obj,scale,T,Pts)
  % obj . plot_dynamics (scale, T, Pts) - Plot spinnet dynamics
  %
  % Plot distances as vectorrs from given points to others using (1 - probability)
  % as distance.
  %
  %   scale  - Plot scaling factor for distances
  %   T      - Time vector (optional, default: inf)
  %   Pts    - Excitations for which to track distance from spins. Single nagative
  %            value means plot for that excitation is in the middle. (optinal,
  %            defaul: 1:N).
  %   obj    - Quantum spin network object
  %


  if ~exist('scale','var')
    scale = 1
  end
  if ~isscalar (scale)
    error ('Illegal scale value');
  end

  if ~exist('T','var') || isempty(T)
    T = inf;
  end
  if ~isvector (T) || min(T) < 0
    error ('Illegal T value');
  end

  cent = false;
  if ~exist('Pts','var')
    Pts = [1:obj.N];
  else
    if sum(abs(Pts) > obj.N || Pts == 0) ~= 0
      error ('Illegal Pts value');
    end
    if length(Pts) == 1 && Pts < 0
      cent = true;
      Pts = -Pts;
    end
  end

  Dm = (1 - obj.prob ());

  M = length(Pts);
  cmap = varycolor(M);

  clf;
  hold on;

  ax=[min([obj.pos(1,:) -.5]) max([obj.pos(1,:) .5]) min([obj.pos(2,:) -.5]) max([obj.pos(2,:) .5])] * 1.2;

  title('Transition Probabilities');

  ang=0:0.01:2*pi;
  xp=scale*cos(ang);
  yp=scale*sin(ang);

  for t = T
    Dt = (1 - obj.prob (t));
    clf;
    axis(ax);
    hold on;
    plot (obj.pos(1,:), obj.pos(2,:), '*b');
    plot (obj.pos(1,Pts), obj.pos(2,Pts), 'or');
    for p = 1:M
      for q = 1:obj.N
        if cent
          if strcmpi(obj.type,'chain')
            s = [obj.pos(1,Pts(p)); .5];
          else
            s = [0; 0];
          end
          v = obj.pos(:,q) - s;
          v = scale * v / norm(v);
          t = s + Dm(Pts(p),q) * v;
          r = s + Dt(Pts(p),q) * v;
          plot([t(1), r(1)], [t(2), r(2)], '-', 'color', cmap(p,:), 'LineWidth', 1);
          plot([s(1), t(1)], [s(2), t(2)], '-', 'color', cmap(p,:)*.8, 'LineWidth', 2);
          plot(s(1)+xp,s(2)+yp, '-.', 'color', cmap(p,:));
        else
          s = obj.pos(:,Pts(p));
          if Pts(p) == q
            plot(s(1)+Dt(Pts(p),q)*xp,s(2)+Dt(Pts(p),q)*yp, ':', 'color', cmap(p,:));
            plot(s(1)+Dm(Pts(p),q)*xp,s(2)+Dm(Pts(p),q)*yp, '-', 'color', cmap(p,:) * .8);
            plot(s(1)+xp,s(2)+yp, '-.', 'color', cmap(p,:));
          else
            if strcmpi(obj.type,'chain')
              if q < Pts(p)
                a = (Pts(p)-q) / (Pts(p)) * pi / 2;
              else
                a = (Pts(p)-q) / (obj.N-Pts(p)) * pi / 2;
              end
              v = [cos(a) -sin(a); sin(a) cos(a)] * [0; scale];
            else
              v = obj.pos(:,q) - obj.pos (:,Pts(p));
              v = scale * v / norm(v);
            end
            t = s + Dm(Pts(p),q) * v;
            r = s + Dt(Pts(p),q) * v;
            plot([t(1), r(1)], [t(2), r(2)], '-', 'color', cmap(p,:), 'LineWidth', 1);
            plot([s(1), t(1)], [s(2), t(2)], '-', 'color', cmap(p,:)*.8, 'LineWidth', 2);
          end
        end
      end
    end
    hold off;
    pause(.01);
  end

end