c      subroutine d01gaf(x,y,n,ans,er,ifail)
      subroutine integran(x,y,n,res,eri,ifail)
c
c     this subroutine integrates a function (y) specified
c     numerically at n points (x), where n is at least 4,
c     over the range x(1) to x(n).  the points need not be
c     equally spaced, but should be distinct and in ascending
c     or descending order.  an error estimate is returned.
c     the method is due to gill and miller.
c
c     nag copyright 1975
c     mark 5 release
c     mark 7 revised ier-154 (dec 1978)
c     mark 11.5(f77) revised. (sept 1985.)
c     .. parameters ..
      character*6       srname
      parameter         (srname='d01gaf')
c     .. scalar arguments ..
      real*4  ans, er
      integer           ifail, n
c     .. array arguments ..
      real*4  x(*), y(*),res(*),eri(*)
c     .. local scalars ..
      real*4  c, d1, d2, d3, h1, h2, h3, h4, r1, r2, r3, r4, s
      integer           i, nn
c     .. external functions ..
      integer           p01abf
      external          p01abf
c     .. executable statements ..
      ans = 0.00
      er = 0.00
c      if (n.ge.4) go to 20
c      ifail = p01abf(ifail,1,srname,0,p01rec)
c      return
c
c     check points are strictly increasing or decreasing
c
   20 h2 = x(2) - x(1)
      do 80 i = 3, n
         h3 = x(i) - x(i-1)
c         if (h2*h3) 40, 60, 80
c   40    ifail = p01abf(ifail,2,srname,0,p01rec)
c         return
c   60    ifail = p01abf(ifail,3,srname,0,p01rec)
c         return
   80 continue
c
c     integrate over initial interval
c
      d3 = (y(2)-y(1))/h2
      h3 = x(3) - x(2)
      d1 = (y(3)-y(2))/h3
      h1 = h2 + h3
      d2 = (d1-d3)/h1
      h4 = x(4) - x(3)
      r1 = (y(4)-y(3))/h4
      r2 = (r1-d1)/(h4+h3)
      h1 = h1 + h4
      r3 = (r2-d2)/h1
      ans = h2*(y(1)+h2*(d3/2.00-h2*(d2/6.00-(h2+2.00*h3)*r3/12.00))
     *      )
      s = -(h2**3)*(h2*(3.00*h2+5.00*h4)+10.00*h3*h1)/60.00
      r4 = 0.00
c
c     integrate over central portion of range
c
      nn = n - 1
      do 120 i = 3, nn
         ans = ans + h3*((y(i)+y(i-1))/2.00-h3*h3*(d2+r2+(h2-h4)*r3)
     *         /12.00)
         c = h3**3*(2.00*h3*h3+5.0*(h3*(h4+h2)+2.00*h4*h2))/120.00
         er = er + (c+s)*r4
         if (i.ne.3) s = c
         if (i.eq.3) s = s + 2.00*c
         eri(i+1)=er-h4**3*r4*(h4*(3.00*h4+5.00*h2)+10.00*h3*(h2+h3+h4))/60.00 + s*r4 
         res(i+1)= ans+h4*(y(i+1)-h4*(r1/2.00+h4*(r2/6.00+(2.00*h3+h4)*r3/12.0 0)))+eri(i+1)
         if (i-n+1) 100, 140, 100
  100    h1 = h2
         h2 = h3
         h3 = h4
         d1 = r1
         d2 = r2
         d3 = r3
         h4 = x(i+2) - x(i+1)
         r1 = (y(i+2)-y(i+1))/h4
         r4 = h4 + h3
         r2 = (r1-d1)/r4
         r4 = r4 + h2
         r3 = (r2-d2)/r4
         r4 = r4 + h1
         r4 = (r3-d3)/r4
  120 continue
c
c     integrate over final interval
c
  140 continue
c      ans = ans + h4*(y(n)-h4*(r1/2.00+h4*(r2/6.00+(2.00*h3+h4)
c     *     *r3/12.00)))
c      er = er - h4**3*r4*(h4*(3.00*h4+5.00*h2)+10.00*h3*(h2+h3+h4))
c     *    /60.00 + s*r4
c		
c      ans = ans + er
      ifail = 0
      return
      end
