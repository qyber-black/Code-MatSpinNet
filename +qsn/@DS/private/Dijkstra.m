function [I,L,N] = Dijkstra(w,s)
  % [I,L,N] = Dijkstra(w,s)
  %
  % This function computes the minimum cost paths from a source node s in
  % a graph defined by its link cost matrix w. This function uses
  % Dijkstra's algorithm.
  %
  % w(i,j) is the cost from node i to node j, w(i,i)=0, and w(i,j)=Inf if
  % there is no direct link from node i to node j.
  %
  % On return, L(v) provides minimum cost from source node s to
  % destination node v. I is the unidirectional incidence matrix of the
  % spanning tree; I(i,j)=1 if there is a direct path from node i to node
  % j in a spanning tree. Otherwise, I(i,j)=0.
  %
  % N is the next node; specifically, N(v) is the next node from source node s
  % in the minimum path from source node s to destination node v.
  %
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019 Edmond Jonckheere, University of Southern California
  % SPDX-License-Identifier: AGPL-3.0-or-later

  % Initialization
  Numnode = size(w,1);
  V = 1:Numnode;
  L = w(s,:);
  I = zeros(Numnode,Numnode);
  N = zeros(1,Numnode);
  N(1,s) = s;
  T = [s];

  % Get Next Node
  for k = 2:Numnode
    VminusT = setdiff(V,T);

    Lj = L(VminusT);
    [Lx,index] = min(Lj);
    vnew = VminusT(index);

    % Find previous node
    for i=1:Numnode
      if ismember(i,T) == 1
        if L(i)+w(i,vnew) == L(vnew)
          vpre = i;
          if i==s
            N(1,vnew) = vnew;
          else
            N(1,vnew) = N(1,vpre);
          end
          break
        end
      end
    end

    % Update least-cost paths
    I(vpre,vnew) = 1;
    T = union(T,vnew);

    for v = VminusT
      if L(v) > L(vnew)+w(vnew,v)
        L(v) = L(vnew)+w(vnew,v);
      end
    end

  end
  
end
