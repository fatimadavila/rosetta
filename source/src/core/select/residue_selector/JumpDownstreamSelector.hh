// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/JumpDownstreamSelector.hh
/// @brief  The JumpDownstreamSelector selects residues downstream of a given jump in a FoldTree
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

#ifndef INCLUDED_core_select_residue_selector_JumpDownstreamSelector_HH
#define INCLUDED_core_select_residue_selector_JumpDownstreamSelector_HH

// Unit headers
#include <core/select/residue_selector/JumpDownstreamSelector.fwd.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/jump_selector/JumpSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

/// @brief The JumpDownstreamSelector returns a ResidueSubset, i.e. a utility::vector1< bool > containing
/// 'true' for residue positions which lie downstream of a given jump in the FoldTree. The jump is
/// specified by an integer index.
class JumpDownstreamSelector : public ResidueSelector {
public:
	// derived from base class
	JumpDownstreamSelector();

	/// @brief Copy constructor
	///
	JumpDownstreamSelector( JumpDownstreamSelector const &src);

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	ResidueSelectorOP clone() const override;

	JumpDownstreamSelector( int jump );
	JumpDownstreamSelector( core::select::jump_selector::JumpSelectorCOP sele );

	~JumpDownstreamSelector() override;

	ResidueSubset apply( core::pose::Pose const & pose ) const override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	) override;

	std::string
	get_name() const override;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	//unit-specific
	/**
	* @brief sets the string by which residues are selected
	*/
	void set_jump( int jump );

	///@brief this overrides the integer setting if this is not null
	void
	set_jump_selector( core::select::jump_selector::JumpSelectorCOP sele ){
		jump_selector_ = sele;
	}

private: // data members
	int jump_ = 0; //this is used iff jump_selector_ == nullptr
	core::select::jump_selector::JumpSelectorCOP jump_selector_ = nullptr;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} //namespace residue_selector
} //namespace select
} //namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_residue_selector_JumpDownstreamSelector )
#endif // SERIALIZATION


#endif
