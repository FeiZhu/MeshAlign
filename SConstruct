#################################################################
# Scons script for MeshAlign Project
# @author: Fei Zhu, 05/28/2015 
#################################################################

#IMPORTS
import fnmatch
import os
import platform
from os.path import basename
from glob import glob

######################## EDIT HERE! #############################
#BUILD TYPE
#build_type='Release'
build_type='Debug'

#BUILD MSVC PROJECTS FOR WINDOWS
build_msvc=True
#build_msvc=False

#PROJECT NAME: SET NAME OF YOUR PROJECT
project_name='MeshAlign'

#USE OPENGL(GL,GLU,GLUT,GLUI) OR NOT
#use_gl=True
use_gl=False

#SRC PATH: SET ROOT DIRECTORY OF YOUR CODE
src_root_path=os.getcwd()

#IGNORED SRC PATH: SET PATH TO DIRECTORIES OF SOURCE CODE WHICH YOU
#                  DON'T WANT TO INCLUDE IN THE BUILD
#EXAMPLE: ignored_src_path=['./code1/','./code2/']
ignored_src_path=[]

#PHYSIKA_PATH: SET PATH TO COMPILED
physika_inc_path=['../Physika/Public_Library/include/']
physika_lib_path=['../Physika/Public_Library/lib/']

#ADDITIONAL PATH: SET ADDITIONAL INCLUDE AND LIB PATH
#EXAMPLE: additional_inc_path=['./path1/','../path2/']
additional_inc_path=[]
additional_lib_path=[]

#PHYSIKA LIBS TO LINK: NOTE NECESSARILLY ALL PHYSIKA LIBS
#NOTE: THE ORDER MATTERS
physika_libs=['Physika_GUI','Physika_Dynamics','Physika_Render','Physika_IO','Physika_Geometry','Physika_Core','LodePNG']

#ADDITIONAL LIBS TO LINK
#EXAMPLE: additional_libs=['lib1','lib2']
additional_libs=[]
##################################################################

#########EDIT ONLY IF YOU'ARE AWARE WHAT YOU'RE DOING#############

#OS TYPE
os_name=platform.system()
os_architecture=platform.architecture()[0]

#COMPILER
compiler=''
if os_name in ('Linux','Darwin') or (os_name=='Windows' and build_msvc==False):
   compiler='g++'
else:
   compiler='msvc'

#PHYSIKA IS DEPENDENT ON STACKWALKER ONLY ON WINDOWS WITH MSVC
if compiler=='msvc':
    physika_libs.append('StackWalker')

#FOR STACK TRACE SUPPORT
if os_name in ('Linux','Darwin'):
    additional_libs.append('dl')
elif compiler=='msvc':
    additional_libs.append('Advapi32')

#SRC FILES
src_files=[]
inc_files=[]
for dir,_,_ in os.walk(src_root_path):
    if dir not in ignored_src_path:
       src_files.extend(glob(os.path.join(dir,'*.cpp')))
       inc_files.extend(glob(os.path.join(dir,'*.h')))

#USE OPENGL OR NOT
if use_gl==True:
    physika_inc_path_str=''.join(physika_inc_path)
    additional_inc_path.append(physika_inc_path_str+'Physika_Dependency/OpenGL/')
    if os_name in ('Linux','Darwin'):
        if os_architecture=='32bit':
            additional_libs.append('glui32')
        else:
            additional_libs.append('glui64')
        additional_libs.append('glut')
        additional_libs.append('X11')
        additional_libs.append('GLU')
        additional_libs.append('GL')
    elif os_name=='Windows':
        additional_libs.append('glui32')
        additional_libs.append('glut')
        additional_libs.append('glu32')
        additional_libs.append('opengl32')

#INCLUDE PATH AND LIB PATH
include_path=physika_inc_path+additional_inc_path
lib_path=physika_lib_path+additional_lib_path
libs=physika_libs+additional_libs

#TARGET NAME
target_name=project_name+build_type

#ENVIRONMENT
ENV={'PATH':os.environ['PATH']}
if compiler=='g++':
   CC='g++'
   CXX='g++'
   tools=['gcc','g++','gnulink']
   LINKFLAGS=['-rdynamic']  #FOR STACK TRACE SUPPORT
   if build_type=='Release':
      CCFLAGS=['-O3','-Wall','-fno-strict-aliasing','-std=gnu++0x','-DNDEBUG']
   else:
      CCFLAGS=['-Wall','-std=gnu++0x','-fno-strict-aliasing','-g']
   env=Environment(CC=CC,CXX=CXX,tools=tools,CCFLAGS=CCFLAGS,LINKFLAGS=LINKFLAGS,CPPPATH=include_path,LIBPATH=lib_path,RPATH=lib_path,LIBS=libs,ENV=ENV)
else:
   if build_type=='Release':
      CCFLAGS=['/Ox','/EHsc','/DNDEBUG','/W3']
   else:
      CCFLAGS=['/Od','/Zi','/EHsc','/W3']
   ENV['TMP']=os.environ['TMP']
   if os_architecture=='32bit':
   	arc='x86'
   else:
	arc='amd64'
   env=Environment(CCFLAGS=CCFLAGS,CPPPATH=include_path,LIBPATH=lib_path,RPATH=lib_path,LIBS=libs,ENV=ENV,MSVS_ARCH=arc,TARGET_ARCH=arc)  
   #DEBUG INFORMATION
   if build_type=='Debug':
      env['LINKFLAGS']=['/DEBUG']
      env['CCPDBFLAGS'] = '/Zi /Fd${TARGET}.pdb'

#BUILD
target=env.Program(target_name,src_files)
#GENERATE MSVC PROJECT OPTIONALLY
if compiler=='msvc':
   sln=env.MSVSProject(target=target_name+env['MSVSPROJECTSUFFIX'],srcs=src_files,incs=inc_files,buildtarget=target,variant=build_type)

#WINDOWS WORKAROUND: COPY DLLS TO EXECUTIVE DIRECTORY
if os_name=='Windows':
    for rpath in lib_path:
        for dll_name in os.listdir(rpath):
            if dll_name.endswith('.dll'):
                Command(dll_name, rpath+dll_name, Copy("$TARGET", "$SOURCE"))

#CUSTOMIZE CLEAN OPTION
if compiler=='msvc':
   sln_delete_files=[build_type+'/','obj/']
   for name in os.listdir('./'):
       if name.endswith('.user') or name.endswith('.pdb') or name.endswith('.suo') or name.endswith('.sdf') or name.endswith('.ilk'):
          sln_delete_files.append(name)
   Clean(sln,sln_delete_files)
