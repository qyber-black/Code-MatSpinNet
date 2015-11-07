function plot_graph(obj, type, scale,violations)
  % obj . plot_graph (type, scale, violations) - Plot distances as graph
  %
  % Requires presence of graphviz software.
  %
  %   type        - Mode of graph construction: 'fdp' or 'neato'
  %   scale       - Graph scale
  %   violations  - Only print triangle violations
  %

  S = obj.D >= obj.eps;
  D = (obj.D * S + 1) * scale;

  filename = tempname ();
  fp = fopen ([filename '.dot'], 'w');
  fprintf (fp, 'graph %s%d {\n  clusterrank=local\n', class(obj.src), obj.N);
  fprintf (fp, '  %d [ color=green, style=filled ];\n', 1:obj.N);

  if violations

    [x,T] = obj.check_triangle_inequality ();

    [A B C] = ind2sub(size(T), find (T < -obj.eps));
    X = sort([A B C],2);

    n = length(A);

    colors = floor(varycolor(n) * 255);

    for I = 1:n
      a = X(I,1);
      b = X(I,2);
      c = X(I,3);
      ok = 1;
      for J = 1:I-1
        if a == X(J,1) && b == X(J,2) && c == X(J,3)
          ok = 0;
          break;
        end
      end
      if ok
        if S(a,b)
          fprintf (fp, '  %d -- %d [ len = %f, color="#%02X%02X%02X" ];\n', a, b, D(a,b), colors(I,1), colors(I,2), colors(I,3));
        else
          fprintf (fp, '  %d -- %d [ len = %f, color=orange ];\n', a, b, D(a,b));
        end
        if S(b,c)
          fprintf (fp, '  %d -- %d [ len = %f, color="#%02X%02X%02X" ];\n', b, c, D(b,c), colors(I,1), colors(I,2), colors(I,3));
        else
          fprintf (fp, '  %d -- %d [ len = %f, color=orange ];\n', b, c, D(b,c));
        end
        if S(c,a)
          fprintf (fp, '  %d -- %d [ len = %f, color="#%02X%02X%02X" ];\n', c, a, D(c,a), colors(I,1), colors(I,2), colors(I,3));
        else
          fprintf (fp, '  %d -- %d [ len = %f, color=orange ];\n', c, a, D(c,a));
        end
      end
    end
    col_eq = 'invis';
    col_neq = 'invis';
  else
    col_eq = 'red';
    col_neq = 'blue';
  end

  for k = 1:obj.N
    for l = k+1:obj.N
      if S(k,l)
        fprintf (fp, '  %d -- %d [ len = %f, color=%s ];\n', k, l, D(k,l), col_neq);
      else
        fprintf (fp, '  %d -- %d [ len = %f, color=%s ];\n', k, l, D(k,l), col_eq);
      end
    end
  end

  fprintf (fp, '}\n');
  fclose(fp);

  system([type ' -T png ' filename '.dot -o ' filename '.png']);
  I = imread ([filename '.png'], 'png');
  clf
  imshow (I);

  delete ([filename '.dot']);
  delete ([filename '.png']);

  pause(0.01);

end