# We need to test this function into oblivion because this set generator is the crux of
# the entire haplotype finding algorithm

# Contrived example: we can classify animals into groups based on their diets or their
# limb structure (stance). Often these things go together (e.g. all ungulates are
# herbavores as far as I know), but somethimes they don't.

# Our function doesn't know anything about zoology, so we create a comparator function
# that says animals can be placed into a set if they all have the same diet, or the same
# stance, but not mixed between the two.

# We then test that animals that we know are grouped are returned as groups, and groups
# that shouldn't exist are not returned.
const STANCES = Dict([
    "🐕" => "digitgrade",
    "🐈" => "digitgrade",
    "🐎" => "unguligrade",
    "🐄" => "unguligrade",
    "🐻" => "plantigrade",
    "🐑" => "unguligrade",
    "👨" => "plantigrade",
    "🐋" => NaN,
    "⛄" => NaN,
])
const DIETS = Dict([
    "🐕" => "carnivore",
    "🐈" => "carnivore",
    "🐎" => "herbavore",
    "🐄" => "herbavore",
    "🐻" => "carnivore",
    "🐑" => "herbavore",
    "👨" => "omnivore",
    "🐋" => "carnivore",
    "⛄" => NaN,
])
const ANIMALS = collect(keys(STANCES))
const CARNIVORES = ["🐕", "🐈", "🐻", "🐋"]
const HERBAVORES = ["🐎", "🐄", "🐑"]
const PLANTIGRADES = ["🐻", "👨"]
const SNOWMAN = ["⛄"]

function set_equals(set)
    return all(p -> STANCES[p] == STANCES[first(set)], set) ||
           all(p -> DIETS[p] == DIETS[first(set)], set)
end #function

function vector_content_equals(a, b)
    return issubset(a, b) && issubset(b, a)
end #function

@testset "findset" begin
    # Test for known groups
    @test any(s -> vector_content_equals(s, CARNIVORES), findset(ANIMALS, set_equals))
    @test any(s -> vector_content_equals(s, HERBAVORES), findset(ANIMALS, set_equals))
    @test any(s -> vector_content_equals(s, PLANTIGRADES), findset(ANIMALS, set_equals))
    @test any(s -> vector_content_equals(s, SNOWMAN), findset(ANIMALS, set_equals))

    # Test against known bad groups
    @test any(s -> !(issubset(["🐋", "⛄"], s)), findset(ANIMALS, set_equals))
    @test any(s -> !(issubset(["🐕", "🐎"], s)), findset(ANIMALS, set_equals))
    @test any(s -> !(issubset(["🐎", "👨"], s)), findset(ANIMALS, set_equals))
    @test any(s -> !(issubset(["🐋", "👨"], s)), findset(ANIMALS, set_equals))

    # Test to make sure nulls did not escape into multiple sets
    @test !any(s -> length(s) > 1, filter(s -> "⛄" in s, findset(ANIMALS, set_equals)))
end #testset
