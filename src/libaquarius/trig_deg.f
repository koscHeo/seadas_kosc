      real*4 function sind(degrees)

      REAL(4), PARAMETER :: PI=3.14159265358979323846D0
      real*4 degrees
      sind = sin(degrees * pi / 180)
      return
      end


      real*8 function dsind(degrees)

      REAL(8), PARAMETER :: PI=3.14159265358979323846D0
      real*8 degrees
      dsind = sin(degrees * pi / 180)
      return
      end


      real*4 function cosd(degrees)

      REAL(4), PARAMETER :: PI=3.14159265358979323846D0
      real*4 degrees
      cosd = cos(degrees * pi / 180)
      return
      end


      real*8 function dcosd(degrees)

      REAL(8), PARAMETER :: PI=3.14159265358979323846D0
      real*8 degrees
      dcosd = cos(degrees * pi / 180)
      return
      end


      real*4 function tand(degrees)

      REAL(4), PARAMETER :: PI=3.14159265358979323846D0
      real*4 degrees
      tand = tan(degrees * pi / 180)
      return
      end


      real*8 function dtand(degrees)

      REAL(8), PARAMETER :: PI=3.14159265358979323846D0
      real*8 degrees
      dtand = tan(degrees * pi / 180)
      return
      end


      real*4 function asind(value)

      REAL(4), PARAMETER :: PI=3.14159265358979323846D0
      real*4 value
      asind = asin(value) * 180 / pi
      return
      end


      real*8 function dasind(value)

      REAL(8), PARAMETER :: PI=3.14159265358979323846D0
      real*8 value
      dasind = asin(value) * 180 / pi
      return
      end


      real*4 function acosd(value)

      REAL(4), PARAMETER :: PI=3.14159265358979323846D0
      real*4 value
      acosd = acos(value) * 180 / pi
      return
      end


      real*8 function dacosd(value)

      REAL(8), PARAMETER :: PI=3.14159265358979323846D0
      real*8 value
      dacosd = acos(value) * 180 / pi
      return
      end


      real*4 function atand(value)

      REAL(4), PARAMETER :: PI=3.14159265358979323846D0
      real*4 value
      atand = atan(value) * 180 / pi
      return
      end


      real*8 function datand(value)

      REAL(8), PARAMETER :: PI=3.14159265358979323846D0
      real*8 value
      datand = atan(value) * 180 / pi
      return
      end


      real*4 function atan2d(value1, value2)

      REAL(4), PARAMETER :: PI=3.14159265358979323846D0
      real*4 value1, value2
      atan2d = atan2(value1, value2) * 180 / pi
      return
      end

      real*8 function datan2d(value1, value2)

      REAL(8), PARAMETER :: PI=3.14159265358979323846D0
      real*8 value1, value2
      datan2d = atan2(value1, value2) * 180 / pi
      return
      end



