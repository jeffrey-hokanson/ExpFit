% function [e,p] = marriage_norm(a,b)
% Given two vectors a and b find the permutation p such that e = norm(a-b(p)) is minimized.
%
% This closely follows the pseudocode on the wikipedia, where a is the man and b is the woman
% (c) Jeffrey Hokanson 1 October 2013, released under GPLv2
function [e,p] = marriage_norm(a,b)
comment = 0;

% ensure that a is always the longest vector
if length(b) < length(a)
	c = b;
	b = a;
	a = c;
	clear c;
end

ma = length(a);
mb = length(b);


% construct preferences for each a
a_pref = zeros(mb,ma);
for j = 1:ma
	[trash,I] = sort(abs(a(j) - b));
	a_pref(:,j) = I;
end

a_engaged = zeros(ma,1);
b_engaged = zeros(mb,1);

while 1
	% find man who hasn't proposed yet
	[t,j] = min(a_engaged);
	if t>0	% if every man has proposed, exit
		p = a_engaged;
		e = norm(a-b(p));
		return
	end
	% if not, now a(j) will propose to b(k)
	% find the most prefered girl j hasn't proposed to
	prefnum = min(find(~isnan(a_pref(:,j))));
	k = a_pref(prefnum,j);

	if comment==1, fprintf('%d proposes to %d ',j,k); end

	if b_engaged(k) == 0 % if the woman isn't engaged, then these two provisionally couple
		b_engaged(k) = j;
		a_engaged(j) = k;
		if comment==1, fprintf('and they become engaged\n'); end
	else % if now, we'll see if the man will be jilted
		fiancee = b_engaged(k);
		if abs(a(fiancee)-b(k))> abs(a(j)-b(k))
			% if the woman likes a(j) more
			b_engaged(k) = j;
			a_engaged(j) = k;
			a_engaged(fiancee) = 0; % ah, the joys of being jilted
			if comment==1, fprintf('and %d replaces %d as her fiancee\n',j,fiancee); end
		else
			if comment==1, fprintf('and woman %d prefers her current fiancee %d\n',k,fiancee); end
			
			% if the woman loves the one she's with, nothing changes
		end
	end
	% in any event, the man has now proposed to this woman, so he can't do so again
	a_pref(prefnum,j) = NaN;
end
