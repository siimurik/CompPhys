!	=============================================================
!	Program to compute roots of a 2nd order polynomial
!	Tasks: Input from user,logical statements,
!	use of functions, stop
!	Accuracy in floating point arithmetic
!	e.g. IF(x.eq.0.0)
!	Tests:		a,	b,	c	=	1		 2		  3					D=	-8.
!				a,	b,	c	=	1		-8		 16					D=	0.
!				a,	b,	c	=	1		-1		 -2					D=	9.
!				a,	b,	c	=	2.3	-2.99		-16.422			D=	160.022507
!				a,	b,	c	=	6.8	-58.48		125.73199		D=	0.000204147349
!				a,	b,	c	=	6.8	-58.48		125.732			D=	-2.10891725E-04
!	=============================================================
program trionymo
implicit none
real(8)	::	a,	b,	c, D
real(8)	::	x1,	x2
real(8)	::	Discriminant

print *, 'Enter a, b, c: '
read *, a, b, c 

!Test if we have a well defined polynomial of 2nd degree:
if (a.eq. 0.0) then 
print *, 'trionymo:		a=0'
print *, 'No longer a 2nd order polynomial'
stop
endif

!Compute the discriminant
D = Discriminant(a, b, c)
print *, 'Discriminant D =', D 

!Compute the roots in each case: D>0, D=0, D<0 (noroots)
if (D.gt.0.0)	then
call roots(a, b, c, x1, x2)
print *, 'Roots:				x1 =', x1, 'x2 =', x2
else if (D.eq.0.0) then
call roots(a, b, c, x1, x2)
print *, 'Double Root:		x1=', x1
else
print *, 'Complex roots'
print*, 'z1 =', -b/(2*a),'+i',abs(sqrt(abs(d))/(2*a)), '&'
print*, 'z2 =', -b/(2*a),'-i',abs(sqrt(abs(d))/(2*a))
endif 

end program trionymo

!	=============================================================
!	This is the function that computes the discriminant
!	A function returns a value. This value is assigned with the
!	statement:
!	Discriminant = <value>
!	i.e. we simply assign anywhere in the program a variable with
!	the name of the function. 
!	=============================================================
real(8) function Discriminant(a, b, c)
implicit none
real(8)		:: a, b ,c

Discriminant = b**2 - 4.0 * a * c

end function Discriminant

!	=============================================================
!	The subroutine that computes the roots.
!	=============================================================
subroutine roots(a, b, c, x1, x2)
implicit none
real(8)		:: a, b, c
real(8)		:: x1, x2
real(8)		:: D, Discriminant

if (a .eq. 0.0) stop  'roots:		a=0'

D = Discriminant(a, b, c)
if (D.ge.0.0) then
D = sqrt(D)
else
print *, 'roots: Sorry, cannot compute roots, D<0 =', D
stop
endif

x1 = (-b + D)/(2.0*a)
x2 = (-b - D)/(2.0*a)

end subroutine roots
