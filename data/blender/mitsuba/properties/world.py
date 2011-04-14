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

from .. import MitsubaAddon

from extensions_framework import declarative_property_group

def MediumParameter(attr, name):
	return [
		{
			'attr': '%s_medium' % attr,
			'type': 'string',
			'name': '%s_medium' % attr,
			'description': '%s medium; blank means vacuum' % attr,
			'save_in_preset': True
		},
		{
			'type': 'prop_search',
			'attr': attr,
			'src': lambda s,c: s.scene.mitsuba_media,
			'src_attr': 'media',
			'trg': lambda s,c: c.mitsuba_material,
			'trg_attr': '%s_medium' % attr,
			'name': name
		}
	]


@MitsubaAddon.addon_register_class
class mitsuba_medium_data(declarative_property_group):
	'''
	Storage class for Mitsuba medium data. The
	mitsuba_media object will store 1 or more of
	these in its CollectionProperty 'media'.
	'''

	ef_attach_to = []	# not attached
	
	controls = [
		'type', 'g', 'densityMultiplier', 'sigmaT', 'albedo'
	]

	properties = [
		{
			'type': 'enum',
			'attr': 'type',
			'name': 'Type',
			'items': [
				('homogeneous', 'Homogeneous', 'homogeneous'),
				('heterogeneous', 'Heterogeneous', 'heterogeneous'),
			],
			'save_in_preset': True
		},
		{
			'type': 'float',
			'attr': 'g',
			'name': 'Asymmetry',
			'description': 'Scattering asymmetry RGB. -1 means back-scattering, 0 is isotropic, 1 is forwards scattering.',
			'default': 0.0,
			'min': -1.0,
			'soft_min': -1.0,
			'max': 1.0,
			'soft_max': 1.0,
			'precision': 4,
			'save_in_preset': True
		},
		{
			'type': 'float',
			'attr': 'densityMultiplier',
			'name': 'Density',
			'description': 'In conjunction with the scattering and absorption coefficients, this number determines the optical density of the medium',
			'default': 1.0,
			'min': 0,
			'max': 10000,
			'precision': 4,
			'save_in_preset': True
		},
		{
			'type': 'float_vector',
			'attr': 'sigmaT',
			'name' : 'Extinction',
			'description' : 'Extinction due to scattering and absorption. Please ' +
				'keep these value roughly equal across color channels (or expect noise).',
			'default' : (1.0, 1.0, 1.0),
			'min': 0.0,
			'max': 1.0,
			'expand' : False,
			'save_in_preset': True
		},
		{
			'type': 'float_vector',
			'attr': 'albedo',
			'subtype': 'COLOR',
			'name' : 'Single-scattering albedo',
			'description' : 'Specifies the albedo of a single scattering interaction',
			'default' : (0.8, 0.8, 0.8),
			'min': 0.0,
			'max': 1.0,
			'expand' : False,
			'save_in_preset': True
		}
	]


@MitsubaAddon.addon_register_class
class mitsuba_media(declarative_property_group):
	'''
	Storage class for Mitsuba Material media.
	'''
	
	ef_attach_to = ['Scene']
	
	controls = [
		'media_select',
		['op_vol_add', 'op_vol_rem']
	]
	
	visibility = {}
	
	properties = [
		{
			'type': 'collection',
			'ptype': mitsuba_medium_data,
			'name': 'media',
			'attr': 'media',
			'items': [
				
			]
		},
		{
			'type': 'int',
			'name': 'media_index',
			'attr': 'media_index',
		},
		{
			'type': 'template_list',
			'name': 'media_select',
			'attr': 'media_select',
			'trg': lambda sc,c: c.mitsuba_media,
			'trg_attr': 'media_index',
			'src': lambda sc,c: c.mitsuba_media,
			'src_attr': 'media',
		},
		{
			'type': 'operator',
			'attr': 'op_vol_add',
			'operator': 'mitsuba.medium_add',
			'text': 'Add',
			'icon': 'ZOOMIN',
		},
		{
			'type': 'operator',
			'attr': 'op_vol_rem',
			'operator': 'mitsuba.medium_remove',
			'text': 'Remove',
			'icon': 'ZOOMOUT',
		},
	]
