#include "RecombChannel.h"

#include <array>

namespace {

constexpr HadronState no_alternative{0, 0.0};

constexpr std::array<RecombChannel, 24> charm_channels{{
    {RecombChannelId::D, ConstituentSystem::Meson, LightFlavor::UpDown, LightFlavor::UpDown,
     OrbitalState::S, {421, 1.86}, {411, 1.87}},
    {RecombChannelId::DStar, ConstituentSystem::Meson, LightFlavor::UpDown, LightFlavor::UpDown,
     OrbitalState::S, {423, 2.01}, {413, 2.01}},
    {RecombChannelId::D1_2420, ConstituentSystem::Meson, LightFlavor::UpDown, LightFlavor::UpDown,
     OrbitalState::P, {10423, 2.42}, no_alternative},
    {RecombChannelId::D0Star_2300, ConstituentSystem::Meson, LightFlavor::UpDown, LightFlavor::UpDown,
     OrbitalState::P, {10411, 2.30}, no_alternative},
    {RecombChannelId::D1_2430, ConstituentSystem::Meson, LightFlavor::UpDown, LightFlavor::UpDown,
     OrbitalState::P, {20413, 2.44}, no_alternative},
    {RecombChannelId::D2Star_2460, ConstituentSystem::Meson, LightFlavor::UpDown, LightFlavor::UpDown,
     OrbitalState::P, {415, 2.46}, no_alternative},
    {RecombChannelId::Ds, ConstituentSystem::Meson, LightFlavor::Strange, LightFlavor::UpDown,
     OrbitalState::S, {431, 1.97}, no_alternative},
    {RecombChannelId::DsStar, ConstituentSystem::Meson, LightFlavor::Strange, LightFlavor::UpDown,
     OrbitalState::S, {433, 2.12}, no_alternative},
    {RecombChannelId::Ds1_2536, ConstituentSystem::Meson, LightFlavor::Strange, LightFlavor::UpDown,
     OrbitalState::P, {10433, 2.54}, no_alternative},
    {RecombChannelId::Ds0Star_2317, ConstituentSystem::Meson, LightFlavor::Strange, LightFlavor::UpDown,
     OrbitalState::P, {10431, 2.32}, no_alternative},
    {RecombChannelId::Ds1_2460, ConstituentSystem::Meson, LightFlavor::Strange, LightFlavor::UpDown,
     OrbitalState::P, {20433, 2.46}, no_alternative},
    {RecombChannelId::Ds2Star_2573, ConstituentSystem::Meson, LightFlavor::Strange, LightFlavor::UpDown,
     OrbitalState::P, {435, 2.57}, no_alternative},
    {RecombChannelId::LambdaC, ConstituentSystem::Baryon, LightFlavor::UpDown, LightFlavor::UpDown,
     OrbitalState::S, {4122, 2.286}, no_alternative},
    {RecombChannelId::LambdaC2860, ConstituentSystem::Baryon, LightFlavor::UpDown, LightFlavor::UpDown,
     OrbitalState::S, {4124, 2.86}, no_alternative},
    {RecombChannelId::LambdaCPWave, ConstituentSystem::Baryon, LightFlavor::UpDown, LightFlavor::UpDown,
     OrbitalState::P, {14122, 2.60}, no_alternative},
    {RecombChannelId::SigmaC, ConstituentSystem::Baryon, LightFlavor::UpDown, LightFlavor::UpDown,
     OrbitalState::S, {4212, 2.46}, no_alternative},
    {RecombChannelId::SigmaCStar, ConstituentSystem::Baryon, LightFlavor::UpDown, LightFlavor::UpDown,
     OrbitalState::S, {4214, 2.52}, no_alternative},
    {RecombChannelId::SigmaCPWave, ConstituentSystem::Baryon, LightFlavor::UpDown, LightFlavor::UpDown,
     OrbitalState::P, {14212, 2.80}, no_alternative},
    {RecombChannelId::XiCGroup1, ConstituentSystem::Unsupported, LightFlavor::UpDown, LightFlavor::Strange,
     OrbitalState::S, no_alternative, no_alternative},
    {RecombChannelId::XiCGroup2, ConstituentSystem::Unsupported, LightFlavor::UpDown, LightFlavor::Strange,
     OrbitalState::S, no_alternative, no_alternative},
    {RecombChannelId::XiCGroup3, ConstituentSystem::Unsupported, LightFlavor::UpDown, LightFlavor::Strange,
     OrbitalState::P, no_alternative, no_alternative},
    {RecombChannelId::OmegaCGroup1, ConstituentSystem::Unsupported, LightFlavor::Strange, LightFlavor::Strange,
     OrbitalState::S, no_alternative, no_alternative},
    {RecombChannelId::OmegaCGroup2, ConstituentSystem::Unsupported, LightFlavor::Strange, LightFlavor::Strange,
     OrbitalState::S, no_alternative, no_alternative},
    {RecombChannelId::OmegaCGroup3, ConstituentSystem::Unsupported, LightFlavor::Strange, LightFlavor::Strange,
     OrbitalState::P, no_alternative, no_alternative},
}};

constexpr std::array<RecombChannel, 6> bottom_channels{{
    {RecombChannelId::B, ConstituentSystem::Meson, LightFlavor::UpDown, LightFlavor::UpDown,
     OrbitalState::S, {521, 5.28}, {511, 5.28}},
    {RecombChannelId::BStar, ConstituentSystem::Meson, LightFlavor::UpDown, LightFlavor::UpDown,
     OrbitalState::S, {513, 5.32}, no_alternative},
    {RecombChannelId::B1, ConstituentSystem::Meson, LightFlavor::UpDown, LightFlavor::UpDown,
     OrbitalState::P, {10513, 5.72}, no_alternative},
    {RecombChannelId::B0Star, ConstituentSystem::Meson, LightFlavor::UpDown, LightFlavor::UpDown,
     OrbitalState::P, {10511, 5.72}, no_alternative},
    {RecombChannelId::B1Star, ConstituentSystem::Meson, LightFlavor::UpDown, LightFlavor::UpDown,
     OrbitalState::P, {20513, 5.75}, no_alternative},
    {RecombChannelId::B2Star, ConstituentSystem::Meson, LightFlavor::UpDown, LightFlavor::UpDown,
     OrbitalState::P, {515, 5.75}, no_alternative},
}};

} // namespace

const RecombChannel* findRecombChannel(
    int HQ_pid, std::size_t channel_index) noexcept
{
  if (HQ_pid == 4 && channel_index < charm_channels.size()) {
    return &charm_channels[channel_index];
  }
  if (HQ_pid == 5 && channel_index < bottom_channels.size()) {
    return &bottom_channels[channel_index];
  }
  return nullptr;
}
