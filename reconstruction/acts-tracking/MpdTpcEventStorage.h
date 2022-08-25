// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

#include <Acts/Utilities/Logger.hpp>

#include <cassert>
#include <memory>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <unordered_map>

namespace Mpd::Tpc {

/// @brief Container to store event data (based on Acts examples).
class EventStorage final {
public:
  EventStorage(Acts::Logging::Level level = Acts::Logging::INFO):
      m_logger(Acts::getDefaultLogger("EventStorage", level)) {}

  EventStorage(const EventStorage&) = delete;
  EventStorage &operator=(const EventStorage&) = delete;

  /// Stores an object.
  template <typename T>
  void add(const std::string &name, T &&object) {
    ACTS_DEBUG("Adding object '" << name << "'");
    assert(!name.empty() && "Object has an empty name");
    assert(m_store.count(name) == 0 && "Object already exists");
    m_store.emplace(name,
                    std::make_unique<HolderT<T>>(std::forward<T>(object)));
    ACTS_VERBOSE("Added object '" << name << "'");
  }

  /// Gets access to a stored object.
  template <typename T>
  const T &get(const std::string &name) const {
    ACTS_DEBUG("Getting object '" << name << "'");
    auto i = m_store.find(name);
    assert(i != m_store.end() && "Object does not exist");

    const IHolder *holder = i->second.get();
    assert(typeid(T) == holder->type());

    ACTS_VERBOSE("Retrieved object '" << name << "'");
    return reinterpret_cast<const HolderT<T>*>(holder)->value;
  }

  /// Checks whether there is an object w/ the given name.
  bool exists(const std::string &name) const {
    return m_store.find(name) != m_store.end();
  }

  /// Clears the storage.
  void clear() {
    m_store.clear();
  }

private:
  // Type-erased value holder for move-constructible types
  struct IHolder {
    virtual ~IHolder() = default;
    virtual const std::type_info &type() const = 0;
  };

  template<typename T, typename =
           std::enable_if_t<std::is_nothrow_move_constructible<T>::value>>
  struct HolderT : public IHolder {
    T value;
    HolderT(T &&v) : value(std::move(v)) {}
    const std::type_info &type() const { return typeid(T); }
  };

  std::unique_ptr<const Acts::Logger> m_logger;
  std::unordered_map<std::string, std::unique_ptr<IHolder>> m_store;

  const Acts::Logger &logger() const { return *m_logger; }
};

} // namespace Mpd::Tpc
