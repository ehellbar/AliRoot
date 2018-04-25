#ifndef ALIEMCALTRIGGERDCSCONFIG_H
#define ALIEMCALTRIGGERDCSCONFIG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//________________________________________________
/// \class AliEMCALTriggerDCSConfig
/// \ingroup EMCALbase
/// \brief Trigger DCS Config
///
/// Add comment
///
/// \author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
//________________________________________________

#include "TObject.h"
#include "TClonesArray.h"
#include <iosfwd>
#include <string>

class AliEMCALTriggerSTUDCSConfig;
class AliEMCALTriggerTRUDCSConfig;

class AliEMCALTriggerDCSConfig : public TObject 
{
  
public:
  
  AliEMCALTriggerDCSConfig();
  virtual ~AliEMCALTriggerDCSConfig();

  /**
   * @brief Equalty operator
   * 
   * Checks if two DCS configs are the same. For equalty all both DCS configurations 
   * and all TRU configurations must match.
   */
  bool operator==(const AliEMCALTriggerDCSConfig &other) const;

  /**
   * @brief Streaming operator for trigger DCS config
   * 
   * Streaming all TRUs and both STUs.
   * 
   * @param stream Stream used for streaming the DCS config object
   * @param config Object to be streamed
   * @return Streaming operator after streaming DCS config 
   */
  friend std::ostream &operator<<(std::ostream &stream, const AliEMCALTriggerDCSConfig &config);
  
	/**
	 * @brief Serialize object to JSON format
	 * 
	 * @return JSON-serialized trigger DCS config object 
	 */
  std::string ToJSON() const;
  
  void                         SetTRUArr(TClonesArray* const ta)             { fTRUArr    = ta; }
  inline void                  SetSTUObj(AliEMCALTriggerSTUDCSConfig* so, Bool_t isDCAL = false);
  
  TClonesArray*                GetTRUArr()                 const             { return fTRUArr;  }
  
  inline AliEMCALTriggerSTUDCSConfig* GetSTUDCSConfig(Bool_t isDCAL = false) const;
  AliEMCALTriggerTRUDCSConfig*        GetTRUDCSConfig(Int_t iTRU) const      { return (AliEMCALTriggerTRUDCSConfig*)fTRUArr->At(iTRU); }

  /**
   * @brief Check whether TRU is enabled
   * 
   * Enabled-status defined via presence of the TRU in the STU region: TRU
   * is enabled if the corresponding bit is set in the STU region
   * 
   * @param iTRU    Global index of the TRU to be checked
   * @return true   TRU is enabled 
   * @return false  TRU is not enabled
   */
  bool IsTRUEnabled(int iTRU) const;

private:
  
  AliEMCALTriggerDCSConfig           (const AliEMCALTriggerDCSConfig &cd); // Not implemented
  AliEMCALTriggerDCSConfig &operator=(const AliEMCALTriggerDCSConfig &cd); // Not implemented
  
  TClonesArray*                fTRUArr;   ///< TRU array
  AliEMCALTriggerSTUDCSConfig* fSTUObj;   ///< STU
  AliEMCALTriggerSTUDCSConfig* fSTUDCAL;  ///< STU of DCAL
  
  /// \cond CLASSIMP
  ClassDef(AliEMCALTriggerDCSConfig,2) ;
  /// \endcond
  
};

void AliEMCALTriggerDCSConfig::SetSTUObj(AliEMCALTriggerSTUDCSConfig* so, Bool_t isDCAL) {
  if(isDCAL) fSTUDCAL = so;
  else       fSTUObj  = so;                                                              }

AliEMCALTriggerSTUDCSConfig* AliEMCALTriggerDCSConfig::GetSTUDCSConfig(Bool_t isDCAL) const {
  if(isDCAL) return fSTUDCAL;
  return fSTUObj;                                                                           }

#endif //ALIEMCALTRIGGERDCSCONFIG_H

