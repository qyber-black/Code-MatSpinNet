function [I,L,N] = BellmanFord(w,s)
  % [I,L,N] = Bellmanford(w,s)
  %
  % This function computes the minimum cost paths from a source node s in
  % a graph defined by its link cost matrix w. This function uses the
  % Bellman-Ford algorithm.
  %
  % w(i,j) is the cost from node i to node j, w(i,i)=0, and w(i,j)=Inf if
  % there is no direct link from node i to node j.
  %
  % L(h,v) is the minimum communication cost from source node s to
  % destination node v, under the constraint that the path from source
  % node s to destination node v has no more than (h-1) links. The last
  % nonvanishing row of L contains the relevant minimum costs.
  %
  % On return, L(v) provides minimum cost from source node s to
  % destination node v. I is the unidirectional incidence matrix of the
  % spanning tree; I(i,j)=1 if there is a direct path from node i to
  % node j in a spanning tree. Otherwise, I(i,j)=0.
  %
  % N is the next node; specifically, N(v) is the next node from source
  % node s in the minimum path from source node s to destination node v.
  %
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019 Edmond Jonckheere, University of Southern California
  % SPDX-License-Identifier: AGPL-3.0-or-later

  % Initialization
  Numnode = size(w,1);
  I = zeros(Numnode,Numnode);
  N = zeros(1,Numnode);
  N(1,s) = s;
  Lh = zeros(Numnode,Numnode);
  Lh(1,:) = Inf;
  Lh(:,s) = 0;

  % Update
  for k = 2:Numnode
    for v = 1:Numnode
      if v ~= s
        [Lh(k,v),vpre] = min(Lh(k-1,:)+(w(:,v))');
        if v ~= vpre
          I(:,v) = 0;
          I(vpre,v) = 1;
        end
        if vpre == s
          N(1,v) = v;
        else
          N(1,v) = N(1,vpre);
        end
      end
    end
    if all((Lh(k,:)-Lh(k-1,:))==0)
      L = Lh(k,:);
      break
    end
  end

end
