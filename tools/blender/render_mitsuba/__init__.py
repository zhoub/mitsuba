# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

bl_addon_info = {
	"name": "Mitsuba",
	"author": "Wenzel Jakob",
	"version": (0, 1),
	"blender": (2, 5, 5),
	"api": 31667,
	"location": "Render > Engine > Mitsuba",
	"description": "Basic Mitsuba integration for Blender",
	"warning": "",
	"wiki_url": "http://wiki.blender.org/index.php/Extensions:2.5/Py/"\
		"Scripts/Render/Mitsuba",
	"tracker_url": "https://www.mitsuba-renderer.org/bugtracker/projects/mitsuba",
	"category": "Render"}


if "bpy" in locals():
	reload(render)
	reload(ui)

else:
	import bpy
	from bpy.props import *
	from render_mitsuba import render
	from render_mitsuba import ui

def register():
	Scene = bpy.types.Scene
	Scene.mts_path = bpy.props.StringProperty(
		name = "mts_path",
		description="Full path to the 'mitsuba' executable",
		default = "", subtype = "FILE_PATH", maxlen=1024)
	Scene.mts_gui = bpy.props.BoolProperty(
		name = "mts_gui",
		description="Should the Mitsuba user interface be launched?",
		default = False)

def unregister():
	del Scene.mts_path

if __name__ == "__main__":
	register()
