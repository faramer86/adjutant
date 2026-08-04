// This library is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.

#include "plink2_zstfile.h"

#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
namespace plink2 {
#endif
// zstd prefix_unknown error string
const char kShortErrZstdPrefixUnknown[] = "Unknown frame descriptor";
const char kShortErrZstdCorruptionDetected[] = "Corrupted block detected";

const char kShortErrZstAlreadyOpen[] = "ZstRfileOpen can't be called on an already-open file";
#ifdef __cplusplus
}
#endif
