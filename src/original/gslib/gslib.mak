# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

!IF "$(CFG)" == ""
CFG=gslib - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to gslib - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "gslib - Win32 Release" && "$(CFG)" != "gslib - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "gslib.mak" CFG="gslib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "gslib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "gslib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
F90=fl32.exe

!IF  "$(CFG)" == "gslib - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
OUTDIR=.\Release
INTDIR=.\Release

ALL : "$(OUTDIR)\gslib.lib"

CLEAN : 
	-@erase ".\Release\gslib.lib"
	-@erase ".\Release\acorni.obj"
	-@erase ".\Release\backtr.obj"
	-@erase ".\Release\beyond.obj"
	-@erase ".\Release\blue.obj"
	-@erase ".\Release\chknam.obj"
	-@erase ".\Release\chktitle.obj"
	-@erase ".\Release\cova3.obj"
	-@erase ".\Release\dlocate.obj"
	-@erase ".\Release\dpowint.obj"
	-@erase ".\Release\dsortem.obj"
	-@erase ".\Release\gauinv.obj"
	-@erase ".\Release\gcum.obj"
	-@erase ".\Release\getindx.obj"
	-@erase ".\Release\getz.obj"
	-@erase ".\Release\green.obj"
	-@erase ".\Release\hexa.obj"
	-@erase ".\Release\ksol.obj"
	-@erase ".\Release\ktsol.obj"
	-@erase ".\Release\locate.obj"
	-@erase ".\Release\nscore.obj"
	-@erase ".\Release\numtext.obj"
	-@erase ".\Release\ordrel.obj"
	-@erase ".\Release\picksupr.obj"
	-@erase ".\Release\powint.obj"
	-@erase ".\Release\psfill.obj"
	-@erase ".\Release\psline.obj"
	-@erase ".\Release\pstext.obj"
	-@erase ".\Release\rand.obj"
	-@erase ".\Release\red.obj"
	-@erase ".\Release\resc.obj"
	-@erase ".\Release\scal.obj"
	-@erase ".\Release\setrot.obj"
	-@erase ".\Release\setsupr.obj"
	-@erase ".\Release\sortem.obj"
	-@erase ".\Release\sqdist.obj"
	-@erase ".\Release\srchsupr.obj"
	-@erase ".\Release\strlen.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Ox /I "Release/" /c /nologo
# ADD F90 /Ox /I "Release/" /c /nologo
F90_PROJ=/Ox /I "Release/" /c /nologo /Fo"Release/" 
F90_OBJS=.\Release/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/gslib.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/gslib.lib" 
LIB32_OBJS= \
	"$(INTDIR)/acorni.obj" \
	"$(INTDIR)/backtr.obj" \
	"$(INTDIR)/beyond.obj" \
	"$(INTDIR)/blue.obj" \
	"$(INTDIR)/chknam.obj" \
	"$(INTDIR)/chktitle.obj" \
	"$(INTDIR)/cova3.obj" \
	"$(INTDIR)/dlocate.obj" \
	"$(INTDIR)/dpowint.obj" \
	"$(INTDIR)/dsortem.obj" \
	"$(INTDIR)/gauinv.obj" \
	"$(INTDIR)/gcum.obj" \
	"$(INTDIR)/getindx.obj" \
	"$(INTDIR)/getz.obj" \
	"$(INTDIR)/green.obj" \
	"$(INTDIR)/hexa.obj" \
	"$(INTDIR)/ksol.obj" \
	"$(INTDIR)/ktsol.obj" \
	"$(INTDIR)/locate.obj" \
	"$(INTDIR)/nscore.obj" \
	"$(INTDIR)/numtext.obj" \
	"$(INTDIR)/ordrel.obj" \
	"$(INTDIR)/picksupr.obj" \
	"$(INTDIR)/powint.obj" \
	"$(INTDIR)/psfill.obj" \
	"$(INTDIR)/psline.obj" \
	"$(INTDIR)/pstext.obj" \
	"$(INTDIR)/rand.obj" \
	"$(INTDIR)/red.obj" \
	"$(INTDIR)/resc.obj" \
	"$(INTDIR)/scal.obj" \
	"$(INTDIR)/setrot.obj" \
	"$(INTDIR)/setsupr.obj" \
	"$(INTDIR)/sortem.obj" \
	"$(INTDIR)/sqdist.obj" \
	"$(INTDIR)/srchsupr.obj" \
	"$(INTDIR)/strlen.obj"

"$(OUTDIR)\gslib.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "gslib - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
OUTDIR=.\Debug
INTDIR=.\Debug

ALL : "$(OUTDIR)\gslib.lib"

CLEAN : 
	-@erase ".\Debug\gslib.lib"
	-@erase ".\Debug\acorni.obj"
	-@erase ".\Debug\backtr.obj"
	-@erase ".\Debug\beyond.obj"
	-@erase ".\Debug\blue.obj"
	-@erase ".\Debug\chknam.obj"
	-@erase ".\Debug\chktitle.obj"
	-@erase ".\Debug\cova3.obj"
	-@erase ".\Debug\dlocate.obj"
	-@erase ".\Debug\dpowint.obj"
	-@erase ".\Debug\dsortem.obj"
	-@erase ".\Debug\gauinv.obj"
	-@erase ".\Debug\gcum.obj"
	-@erase ".\Debug\getindx.obj"
	-@erase ".\Debug\getz.obj"
	-@erase ".\Debug\green.obj"
	-@erase ".\Debug\hexa.obj"
	-@erase ".\Debug\ksol.obj"
	-@erase ".\Debug\ktsol.obj"
	-@erase ".\Debug\locate.obj"
	-@erase ".\Debug\nscore.obj"
	-@erase ".\Debug\numtext.obj"
	-@erase ".\Debug\ordrel.obj"
	-@erase ".\Debug\picksupr.obj"
	-@erase ".\Debug\powint.obj"
	-@erase ".\Debug\psfill.obj"
	-@erase ".\Debug\psline.obj"
	-@erase ".\Debug\pstext.obj"
	-@erase ".\Debug\rand.obj"
	-@erase ".\Debug\red.obj"
	-@erase ".\Debug\resc.obj"
	-@erase ".\Debug\scal.obj"
	-@erase ".\Debug\setrot.obj"
	-@erase ".\Debug\setsupr.obj"
	-@erase ".\Debug\sortem.obj"
	-@erase ".\Debug\sqdist.obj"
	-@erase ".\Debug\srchsupr.obj"
	-@erase ".\Debug\strlen.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Z7 /I "Debug/" /c /nologo
# ADD F90 /Z7 /I "Debug/" /c /nologo
F90_PROJ=/Z7 /I "Debug/" /c /nologo /Fo"Debug/" 
F90_OBJS=.\Debug/
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/gslib.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/gslib.lib" 
LIB32_OBJS= \
	"$(INTDIR)/acorni.obj" \
	"$(INTDIR)/backtr.obj" \
	"$(INTDIR)/beyond.obj" \
	"$(INTDIR)/blue.obj" \
	"$(INTDIR)/chknam.obj" \
	"$(INTDIR)/chktitle.obj" \
	"$(INTDIR)/cova3.obj" \
	"$(INTDIR)/dlocate.obj" \
	"$(INTDIR)/dpowint.obj" \
	"$(INTDIR)/dsortem.obj" \
	"$(INTDIR)/gauinv.obj" \
	"$(INTDIR)/gcum.obj" \
	"$(INTDIR)/getindx.obj" \
	"$(INTDIR)/getz.obj" \
	"$(INTDIR)/green.obj" \
	"$(INTDIR)/hexa.obj" \
	"$(INTDIR)/ksol.obj" \
	"$(INTDIR)/ktsol.obj" \
	"$(INTDIR)/locate.obj" \
	"$(INTDIR)/nscore.obj" \
	"$(INTDIR)/numtext.obj" \
	"$(INTDIR)/ordrel.obj" \
	"$(INTDIR)/picksupr.obj" \
	"$(INTDIR)/powint.obj" \
	"$(INTDIR)/psfill.obj" \
	"$(INTDIR)/psline.obj" \
	"$(INTDIR)/pstext.obj" \
	"$(INTDIR)/rand.obj" \
	"$(INTDIR)/red.obj" \
	"$(INTDIR)/resc.obj" \
	"$(INTDIR)/scal.obj" \
	"$(INTDIR)/setrot.obj" \
	"$(INTDIR)/setsupr.obj" \
	"$(INTDIR)/sortem.obj" \
	"$(INTDIR)/sqdist.obj" \
	"$(INTDIR)/srchsupr.obj" \
	"$(INTDIR)/strlen.obj"

"$(OUTDIR)\gslib.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

################################################################################
# Begin Target

# Name "gslib - Win32 Release"
# Name "gslib - Win32 Debug"

!IF  "$(CFG)" == "gslib - Win32 Release"

!ELSEIF  "$(CFG)" == "gslib - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\acorni.for

"$(INTDIR)\acorni.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\backtr.for

"$(INTDIR)\backtr.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\beyond.for

"$(INTDIR)\beyond.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\blue.for

"$(INTDIR)\blue.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\chknam.for

"$(INTDIR)\chknam.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\chktitle.for

"$(INTDIR)\chktitle.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\cova3.for

"$(INTDIR)\cova3.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\dlocate.for

"$(INTDIR)\dlocate.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\dpowint.for

"$(INTDIR)\dpowint.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\dsortem.for

"$(INTDIR)\dsortem.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\gauinv.for

"$(INTDIR)\gauinv.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\gcum.for

"$(INTDIR)\gcum.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\getindx.for

"$(INTDIR)\getindx.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\getz.for

"$(INTDIR)\getz.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\green.for

"$(INTDIR)\green.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\hexa.for

"$(INTDIR)\hexa.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\ksol.for

"$(INTDIR)\ksol.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\ktsol.for

"$(INTDIR)\ktsol.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\locate.for

"$(INTDIR)\locate.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\nscore.for

"$(INTDIR)\nscore.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\numtext.for

"$(INTDIR)\numtext.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\ordrel.for

"$(INTDIR)\ordrel.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\picksupr.for

"$(INTDIR)\picksupr.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\powint.for

"$(INTDIR)\powint.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\psfill.for

"$(INTDIR)\psfill.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\psline.for

"$(INTDIR)\psline.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\pstext.for

"$(INTDIR)\pstext.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\rand.for

"$(INTDIR)\rand.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\red.for

"$(INTDIR)\red.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\resc.for

"$(INTDIR)\resc.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\scal.for

"$(INTDIR)\scal.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\setrot.for

"$(INTDIR)\setrot.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\setsupr.for

"$(INTDIR)\setsupr.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\sortem.for

"$(INTDIR)\sortem.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\sqdist.for

"$(INTDIR)\sqdist.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\srchsupr.for

"$(INTDIR)\srchsupr.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\strlen.for

"$(INTDIR)\strlen.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
# End Target
# End Project
################################################################################
