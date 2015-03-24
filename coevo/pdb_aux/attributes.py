#!/usr/bin/env python
''' attributes -- utilities for making chimera attributes

'''

def mangle_attr_name(attr_name):
    ''' Chimera attributes can't start with digit, underscore, or capital

    '''

    return 'a_' + attr_name

def make_chimera_attributes(attr_df, chain_id):
    ''' Returns string of attribute control and assignment lines

        attr_df: pandas.Data.Frame in wide format, indexed by resnum
        chain_id: chain id of residues involved

    '''

    attr_str = ''

    for attr_name in attr_df:
        attr_str += 'attribute: %s\n' % mangle_attr_name(attr_name)
        attr_str += 'match mode: 1-to-1\n'
        attr_str += 'recipient: residues\n'
        for (resn, attr_val) in attr_df[attr_name].iteritems():
            atom_spec = ':%d.%s' % (resn, chain_id)
            attr_str += '\t%s\t%s\n' % (atom_spec, attr_val)

    return attr_str

