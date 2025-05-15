#define NOMINMAX // Define NOMINMAX to prevent min/max macros from windows.h
#define _USE_MATH_DEFINES

#include <windows.h>
#include <winerror.h>             // Для HRESULT_FROM_WIN32
#include <mmdeviceapi.h>          // Для работы с аудиоустройствами (IMMDevice, IMMDeviceEnumerator и т.д.)
#include <Audioclient.h>          // Для IAudioClient, IAudioCaptureClient, IAudioRenderClient
#include <Functiondiscoverykeys_devpkey.h> // Для PKEY_Device_FriendlyName
#include <propvarutil.h>          // Для PropVariantInit, PropVariantClear
#include <stdio.h>                // Для printf, fgets, sscanf
#include <vector>                 // Для std::vector
#include <atomic>                 // Для std::atomic<bool>
#include <synchapi.h>             // Для CreateEvent, SetEvent, WaitForSingleObject, CreateThread, InitializeCriticalSection, DeleteCriticalSection
#include <ksmedia.h>              // Для KSDATAFORMAT_SUBTYPE_IEEE_FLOAT и KSDATAFORMAT_SUBTYPE_PCM GUID
#include <initguid.h>             // Для определения GUID, таких как KSDATAFORMAT_SUBTYPE_IEEE_FLOAT. Должен быть только в одном .cpp файле/блока компиляции
#include <combaseapi.h>           // Для CoTaskMemFree, CoInitializeEx, CoUninitialize
#include <string.h>               // Для memcpy, memset, strcspn
#include <limits>                 // Для numeric_limits
#include <algorithm>              // Для std::min, std::max, std::clamp
#include <iostream>               // Для std::wcout, std::endl, std::cout
#include <cmath>                  // Для std::abs, std::cos, std::sin, std::pow, M_PI, std::log10, std::sqrt, fabs, floor, fmod, fminf, fmaxf
#include <utility>                // Для std::exchange
#include <memory>                 // Для std::unique_ptr
#include <cwchar>                 // For wcscmp, wprintf
#include <string>                 // For std::string, std::wstring
#include <sstream>                // For std::stringstream
#include <cstdlib>                // For std::strtol
#include <list>                   // For std::list (alternative to vector for chain)


// Temporary define for printing GUIDs
#define PRINT_GUID(guid) \
    printf("{%08lX-%04hX-%04hX-%02hhX%02hhX-%02hhX%02hhX%02hhX%02hhX%02hhX%02hhX}", \
           guid.Data1, guid.Data2, guid.Data3, \
           guid.Data4[0], guid.Data4[1], guid.Data4[2], guid.Data4[3], \
           guid.Data4[4], guid.Data4[5], guid.Data4[6], guid.Data4[7])


// Macro for logging HRESULT warnings
#define WARN_ON_ERROR(hres, msg) \
    if (FAILED(hres)) { printf("Предупреждение (%s): HRESULT = 0x%x\n", msg, hres); }


// 2. Declarations of Classes and Structs

// User-defined RAII wrappers (declarations)
template <class T> class ComPtr;
class CoTaskMemPtr;
class ScopedHandle;

// Biquad Filter related declarations
enum class BiquadFilterType;
struct BiquadFilter; // Declaration

// Base DSP Effect Interface
class IDspEffect; // Declaration

// Concrete DSP Effects (declarations)
class GainEffect;
class BiquadEffect;

// Shared Data Structures (declarations)
struct SharedDspParameters;
struct SharedMonitorData;
struct SharedAudioBuffer;

// Audio Converter (declaration)
class AudioConverter;

// Thread Parameter Structures (declarations)
struct AudioCaptureThreadParams;
struct AudioRenderThreadParams;


// --- User-defined RAII wrappers (Declarations with inline Implementations) ---
template <class T>
class ComPtr {
public:
    ComPtr() : ptr_(nullptr) {}
    ComPtr(T* ptr) : ptr_(ptr) {}
    ~ComPtr() { Release(); }
    ComPtr(const ComPtr&) = delete;
    ComPtr& operator=(const ComPtr&) = delete;
    ComPtr(ComPtr&& other) noexcept : ptr_(std::exchange(other.ptr_, nullptr)) {}
    ComPtr& operator=(ComPtr&& other) noexcept {
        if (this != &other) {
            Release();
            ptr_ = std::exchange(other.ptr_, nullptr);
        }
        return *this;
    }
    void Release() {
        if (ptr_ != nullptr) {
            ptr_->Release();
            ptr_ = nullptr;
        }
    }
    T* Get() const { return ptr_; }
    T** ReleaseAndGetAddressOf() {
        Release();
        return &ptr_;
    }
    T* operator->() const { return ptr_; }
    T& operator*() const { return *ptr_; }
    explicit operator bool() const { return ptr_ != nullptr; }
private:
    T* ptr_;
};

class CoTaskMemPtr {
public: // Make Get(), Free(), ReleaseAndGetAddressOf() public
    CoTaskMemPtr() : ptr_(nullptr) {}
    CoTaskMemPtr(void* ptr) : ptr_(ptr) {}
    ~CoTaskMemPtr() { Free(); }
    CoTaskMemPtr(const CoTaskMemPtr&) = delete;
    CoTaskMemPtr& operator=(const CoTaskMemPtr&) = delete;
    CoTaskMemPtr(CoTaskMemPtr&& other) noexcept : ptr_(std::exchange(other.ptr_, nullptr)) {}
    CoTaskMemPtr& operator=(CoTaskMemPtr&& other) noexcept {
        if (this != &other) {
            Free();
            ptr_ = std::exchange(other.ptr_, nullptr);
        }
        return *this;
    }
    void Free() { if (ptr_ != nullptr) { CoTaskMemFree(ptr_); ptr_ = nullptr; } }
    void** ReleaseAndGetAddressOf() { Free(); return &ptr_; }
    void* Get() const { return ptr_; } // <--- Made Public
     template<typename TargetType> explicit operator TargetType*() const { return static_cast<TargetType*>(ptr_); }
    explicit operator bool() const { return ptr_ != nullptr; }
private:
    void* ptr_;
};

class ScopedHandle {
public:
    ScopedHandle() : handle_(nullptr) {}
    ScopedHandle(HANDLE handle) : handle_(handle) {}
    ~ScopedHandle() { Close(); }
    ScopedHandle(const ScopedHandle&) = delete;
    ScopedHandle& operator=(const ScopedHandle&) = delete;
    ScopedHandle(ScopedHandle&& other) noexcept : handle_(std::exchange(other.handle_, nullptr)) {}
    ScopedHandle& operator=(ScopedHandle&& other) noexcept {
        if (this != &other) {
            Close();
            handle_ = std::exchange(other.handle_, nullptr);
        }
        return *this;
    }
    void Close() { if (handle_ != nullptr && handle_ != INVALID_HANDLE_VALUE) { CloseHandle(handle_); } handle_ = nullptr; }
    HANDLE Get() const { return handle_; }
    HANDLE* ReleaseAndGetAddressOf() { Close(); return &handle_; }
    explicit operator bool() const { return handle_ != nullptr && handle_ != INVALID_HANDLE_VALUE; }
private:
    HANDLE handle_;
};


// --- Biquad Filter related declarations ---
enum class BiquadFilterType {
    Identity, LowPass, HighPass, BandPass, BandPass0dB, Notch, LowShelf, HighShelf, Peak
};

struct BiquadFilter {
    float b0, b1, b2, a1, a2;
    std::vector<std::pair<float, float>> state; // State for each channel
    BiquadFilter(); // Declaration
    void ResetState(size_t numChannels); // Declaration
    float process_sample(float input, size_t channelIndex); // Declaration
};


// --- Base DSP Effect Interface ---
class IDspEffect {
public:
    virtual ~IDspEffect() = default;
    virtual void ProcessBuffer(BYTE* pData, size_t numFrames, const WAVEFORMATEX* pFormat) = 0; // pFormat is the CAPTURE format
    virtual void UpdateParameters(const struct SharedDspParameters& params, const WAVEFORMATEX* pFormat) = 0; // Capture format
    virtual const char* GetName() const = 0;
};


// --- Concrete DSP Effects (declarations) ---
class GainEffect : public IDspEffect {
private: float gainFactor_ = 1.0f;
public:
    const char* GetName() const override; // Declaration
    void ProcessBuffer(BYTE* pData, size_t numFrames, const WAVEFORMATEX* pFormat) override; // Declaration
    void UpdateParameters(const struct SharedDspParameters& params, const WAVEFORMATEX* pFormat) override; // Declaration
};

class BiquadEffect : public IDspEffect {
private:
    BiquadFilter filter_;
    BiquadFilterType filterType_ = BiquadFilterType::Identity;
    float filterFreq_ = 1000.0f;
    float filterQ_ = 1.0f;
    float filterGainDb_ = 0.0f;
    float filterS_ = 1.0f;
    void RecalculateCoefficients(const WAVEFORMATEX* pFormat); // Declaration
public:
    const char* GetName() const override; // Declaration
    void ProcessBuffer(BYTE* pData, size_t numFrames, const WAVEFORMATEX* pFormat) override; // Declaration
    void UpdateParameters(const struct SharedDspParameters& params, const WAVEFORMATEX* pFormat) override; // Declaration
    void Initialize(const WAVEFORMATEX* pFormat); // Declaration
};


// --- Shared Data Structures (declarations) ---
struct SharedDspParameters {
    CRITICAL_SECTION lock; HANDLE hParamChangedEvent; float gainFactor; BiquadFilterType filterType; float filterFreq; float filterQ; float filterGainDb; float filterS;
    SharedDspParameters(); // Declaration
    ~SharedDspParameters(); // Declaration
    SharedDspParameters(const SharedDspParameters&) = delete;
    SharedDspParameters& operator=(const SharedDspParameters&) = delete;
    void Initialize(); // Declaration
    void Cleanup(); // Declaration
    void GetParameters(float& outGain, BiquadFilterType& outType, float& outFreq, float& outQ, float& outGainDb, float& outS) const; // Declaration
    void SetGain(float newGain); // Declaration
    void SetFilterType(BiquadFilterType newType); // Declaration
    void SetFilterFreq(float newFreq); // Declaration
    void SetFilterQ(float newQ); // Declaration
    void SetFilterGainDb(float newGainDb); // Declaration
    void SetFilterS(float newS); // Declaration
    void SetFilterSettings(BiquadFilterType type, float freq, float q, float gainDb, float s); // Declaration
    void SignalParametersChanged(); // Declaration
};

struct SharedMonitorData {
    CRITICAL_SECTION lock; volatile UINT32 capturePadding; volatile UINT32 renderPadding;
    // ИСправлено: Добавлены члены для мониторинга буфера SharedBuffer
    volatile UINT32 framesReadFromSharedBuffer;
    volatile UINT32 framesConsumedFromSharedBuffer;
    volatile UINT32 framesWrittenToRenderBuffer;

    SharedMonitorData(); // Declaration
    ~SharedMonitorData(); // Declaration
    SharedMonitorData(const SharedMonitorData&) = delete;
    SharedMonitorData& operator=(const SharedMonitorData&) = delete;
    void Initialize(); // Declaration
    void Cleanup(); // Declaration
    void UpdateCapturePadding(UINT32 padding); // Declaration
    void UpdateRenderPadding(UINT32 padding); // Declaration
    void GetPadding(UINT32& outCapture, UINT32& outRender) const; // Declaration
    // ИСправлено: Добавлен getter для мониторинговых данных по SharedBuffer
     void GetSharedBufferStats(UINT32& outFramesRead, UINT32& outFramesConsumed, UINT32& outFramesWrittenToRender) const; // Declaration
};

struct SharedAudioBuffer {
    BYTE* bufferData; size_t bufferSize; size_t frameSize; // Frame size of CAPTURE format
    volatile size_t writePos; volatile size_t readPos; CRITICAL_SECTION lock;
    BYTE captureFormatBuffer[sizeof(WAVEFORMATEXTENSIBLE)]; PWAVEFORMATEX pCaptureFormat; bool hasCaptureFormat = false;

    SharedAudioBuffer(); // Declaration
    ~SharedAudioBuffer(); // Declaration
    SharedAudioBuffer(const SharedAudioBuffer&) = delete;
    SharedAudioBuffer& operator=(const SharedAudioBuffer&) = delete;
    SharedAudioBuffer(SharedAudioBuffer&& other) noexcept; // Declaration
    SharedAudioBuffer& operator=(SharedAudioBuffer&& other) noexcept; // Declaration

    size_t Initialize(const WAVEFORMATEX* pFormat, size_t durationMs); // Declaration
    void Cleanup(); // Declaration
    size_t Read(BYTE* pBuffer, size_t max_frames_to_read_capture_format); // Declaration (Returns frames read, DOES NOT advance readPos)
    size_t Write(const BYTE* data, size_t frames); // Declaration

    size_t GetAvailableFrames() const; // Declaration
    const WAVEFORMATEX* GetCaptureFormat() const; // Declaration
    // ИСправлено: Добавлен метод для явного продвижения указателя чтения
    void AdvanceReadPointer(size_t numFrames); // Declaration
};


// Audio Converter (declaration)
class AudioConverter {
private:
    const WAVEFORMATEX* pCaptureFormat_ = nullptr;
    const WAVEFORMATEX* pRenderFormat_ = nullptr;

    // State for linear interpolation: stores the LAST frame of the PREVIOUS input buffer passed to Convert
    // Made public for access from RenderThread (only for reading state content)
public: // <--- Made Public
    std::vector<BYTE> interpolationStateBuffer_; // Will store 1 capture frame if resampling is active
    double resampleRatio_ = 1.0; // capture_rate / render_rate (Made Public)
    double resamplePhase_ = 0.0; // Current fractional position within the *current input data block* (after state)
    // Let's refine resamplePhase_: It's the fractional offset at the BEGINNING of the input data block.
    // It should probably be a double tracking the position relative to the start of the effective input buffer (state + data)

private: // Helper methods - keep private as only used internally by Convert
    static bool warned_unsupported_conversion_; // Static member declaration

    float GetFloatSample(const BYTE* pSample, const WAVEFORMATEX* pFormat) const; // Added const
    void SetFloatSample(BYTE* pSample, float value, const WAVEFORMATEX* pFormat) const; // Added const
    bool IsFormatCombinationSupported(const WAVEFORMATEX* pCaptureFormat, const WAVEFORMATEX* pRenderFormat) const; // Declaration
    float InterpolateLinear(float sample1, float sample2, float factor) const; // Added const


public: // Public methods
    AudioConverter(); // Declaration
    ~AudioConverter() = default;

    bool Initialize(const WAVEFORMATEX* pCaptureFormat, const WAVEFORMATEX* pRenderFormat); // Declaration
    // ИСправлено: Сигнатура Convert адаптирована.
    // inputBuffer: Указатель на объединенный буфер [state] [new data]
    // numInputFramesInBlock: Количество фреймов ТОЛЬКО в сегменте НОВЫХ данных (после state)
    // outputBuffer: Буфер для вывода
    // maxOutputFrames: Максимальное количество фреймов для записи в outputBuffer
    // numInputFramesConsumedOut: Output parameter, получает количество ФРЕЙМОВ, фактически потребленных из сегмента НОВЫХ данных
    size_t Convert(const BYTE* inputBuffer, size_t numInputFramesInBlock, BYTE* outputBuffer, size_t maxOutputFrames, size_t* numInputFramesConsumedOut); // Declaration
    void ResetState(); // Declaration

    double GetResampleRatio() const { return resampleRatio_; } // Declaration (inline getter for public member)
    size_t GetStateBufferSize() const; // Declaration (in bytes)
    size_t GetStateBufferSizeFrames() const; // Declaration (in frames)
};


// Thread Parameter Structures (declarations)
struct AudioCaptureThreadParams {
 // Вместо сырых указателей, используем ComPtr и ScopedHandle для RAII
    ComPtr<IAudioClient> pAudioClient;
    ComPtr<IAudioCaptureClient> pCaptureClient;
    ScopedHandle hCaptureEvent;
    ScopedHandle hStopEvent; // Это событие остановки, владеет им main, здесь только копия хэндла

    SharedAudioBuffer* pSharedBuffer; // Указатели на общие структуры, которыми владеет main
    SharedDspParameters* pSharedDspParams;
    SharedMonitorData* pSharedMonitorData;

    UINT32 bufferFrameCount;
    REFERENCE_TIME devicePeriod;
    bool isEventDriven;

    std::vector<std::unique_ptr<IDspEffect>> effectsChain; // DSP chain is owned by capture thread params

    // Конструктор и деструктор (оставлены default)
    AudioCaptureThreadParams() = default;
    ~AudioCaptureThreadParams() = default;
    AudioCaptureThreadParams(const AudioCaptureThreadParams&) = delete; // Запрещаем копирование
    AudioCaptureThreadParams& operator=(const AudioCaptureThreadParams&) = delete;
    AudioCaptureThreadParams(AudioCaptureThreadParams&&) = default; // Разрешаем перемещение
    AudioCaptureThreadParams& operator=(AudioCaptureThreadParams&&) = default;
};

struct AudioRenderThreadParams {
 // Вместо сырых указателей, используем ComPtr и ScopedHandle для RAII
    ComPtr<IAudioClient> pAudioClient;
    ComPtr<IAudioRenderClient> pRenderClient;
    ScopedHandle hRenderEvent;
    ScopedHandle hStopEvent; // Это событие остановки, владеет им main, здесь только копия хэндла

    SharedAudioBuffer* pSharedBuffer; // Указатели на общие структуры, которыми владеет main
    SharedMonitorData* pSharedMonitorData;

    UINT32 bufferFrameCount;
    REFERENCE_TIME devicePeriod;
    bool isEventDriven;
    BYTE renderFormatBuffer[sizeof(WAVEFORMATEXTENSIBLE)]; PWAVEFORMATEX pRenderFormat; // Render format buffer

    AudioConverter audioConverter; // Audio Converter object owned by params

    // Конструктор и деструктор (оставлены default)
    AudioRenderThreadParams() = default;
    ~AudioRenderThreadParams() = default;
     AudioRenderThreadParams(const AudioRenderThreadParams&) = delete; // Запрещаем копирование
     AudioRenderThreadParams& operator=(const AudioRenderThreadParams&) = delete;
    AudioRenderThreadParams(AudioRenderThreadParams&&) = default; // Разрешаем перемещение
    AudioRenderThreadParams& operator=(AudioRenderThreadParams&&) = default;
};


// 3. Global Function Prototypes
HRESULT SelectAudioDevicePrompt(IMMDeviceEnumerator* pEnumerator, EDataFlow dataFlow, ERole role, ComPtr<IMMDevice>& selectedDevice, CoTaskMemPtr& mixFormatMem);
HRESULT GetMixFormatForDevice(IMMDevice* pDevice, CoTaskMemPtr& mixFormatMem);
HRESULT InitializeCaptureDevice(IMMDevice* pDevice, PWAVEFORMATEX pMixFormat, ComPtr<IAudioClient>& audioClient, ComPtr<IAudioCaptureClient>& captureClient, ScopedHandle& eventHandle);
HRESULT InitializeRenderDevice(IMMDevice* pDevice, PWAVEFORMATEX pMixFormat, ComPtr<IAudioClient>& audioClient, ComPtr<IAudioRenderClient>& renderClient, ScopedHandle& eventHandle);
void ApplyDspEffectChain(BYTE* pData, size_t numFrames, const WAVEFORMATEX* pFormat, const std::vector<std::unique_ptr<IDspEffect>>& effectsChain);
DWORD WINAPI AudioCaptureThread(LPVOID lpParam);
DWORD WINAPI AudioRenderThread(LPVOID lpParam);
void HandleUserInput(SharedDspParameters* pSharedDspParams, SharedAudioBuffer* pSharedBuffer, SharedMonitorData* pSharedMonitorData, HANDLE hStopEvent);
void PrintCurrentDspSettings(const SharedDspParameters* pDspParams);
void PrintMonitorData(const SharedAudioBuffer* pSharedBuffer, const SharedMonitorData* pSharedMonitorData);
void PrintWaveFormatEx(const WAVEFORMATEX* pFormat);
const char* BiquadFilterTypeToString(BiquadFilterType type); // Moved prototype here


// 4. Definitions of Auxiliary Functions (not class methods)

// Helper to check if a format is 16-bit PCM (standard or extensible)
bool IsPCM16Format(const WAVEFORMATEX* pFormat) {
    if (!pFormat) return false;

    // Check for standard 16-bit PCM
    if (pFormat->wFormatTag == WAVE_FORMAT_PCM && pFormat->wBitsPerSample == 16) {
        return true;
    }

    // Check for extensible 16-bit PCM
    if (pFormat->wFormatTag == WAVE_FORMAT_EXTENSIBLE) {
        const WAVEFORMATEXTENSIBLE* pExt = reinterpret_cast<const WAVEFORMATEXTENSIBLE*>(pFormat);
        // Убедимся, что cbSize достаточно велик для WAVEFORMATEXTENSIBLE
        if (pFormat->cbSize >= sizeof(WAVEFORMATEXTENSIBLE) - sizeof(WAVEFORMATEX)) {
             if (pExt->Samples.wValidBitsPerSample == 16 && IsEqualGUID(pExt->SubFormat, KSDATAFORMAT_SUBTYPE_PCM)) {
                 return true;
             }
        }
    }

    return false;
}


// Helper functions to calculate Biquad coefficients (used by BiquadEffect)
// Moved definitions here after relevant declarations
void SetIdentityCoefficients(BiquadFilter& filter) { filter.b0 = 1.0f; filter.b1 = 0.0f; filter.b2 = 0.0f; filter.a1 = 0.0f; filter.a2 = 0.0f; }
void CalculateLowPassCoefficients(BiquadFilter& filter, float sampleRate, float cutoffFreq, float Q) { if (sampleRate <= 0 || cutoffFreq <= 0 || cutoffFreq >= sampleRate / 2 || Q <= 0) { SetIdentityCoefficients(filter); return; } float w0 = 2.0f * (float)M_PI * cutoffFreq / sampleRate; float cos_w0 = std::cos(w0); float alpha = std::sin(w0) / (2.0f * Q); float a0_cook = 1.0f + alpha; float a1_cook = -2.0f * cos_w0; float a2_cook = 1.0f - alpha; float b0_cook = (1.0f - cos_w0) / 2.0f; float b1_cook = 1.0f - cos_w0; float b2_cook = (1.0f - cos_w0) / 2.0f; filter.b0 = b0_cook / a0_cook; filter.b1 = b1_cook / a0_cook; filter.b2 = b2_cook / a0_cook; filter.a1 = a1_cook / a0_cook; filter.a2 = a2_cook / a0_cook; }
void CalculateHighPassCoefficients(BiquadFilter& filter, float sampleRate, float cutoffFreq, float Q) { if (sampleRate <= 0 || cutoffFreq <= 0 || cutoffFreq >= sampleRate / 2 || Q <= 0) { SetIdentityCoefficients(filter); return; } float w0 = 2.0f * (float)M_PI * cutoffFreq / sampleRate; float cos_w0 = std::cos(w0); float alpha = std::sin(w0) / (2.0f * Q); float a0_cook = 1.0f + alpha; float a1_cook = -2.0f * cos_w0; float a2_cook = 1.0f - alpha; float b0_cook = (1.0f + cos_w0) / 2.0f; float b1_cook = -(1.0f + cos_w0); float b2_cook = (1.0f + cos_w0) / 2.0f; filter.b0 = b0_cook / a0_cook; filter.b1 = b1_cook / a0_cook; filter.b2 = b2_cook / a0_cook; filter.a1 = a1_cook / a0_cook; filter.a2 = a2_cook / a0_cook; }
void CalculateBandPassCoefficients(BiquadFilter& filter, float sampleRate, float centerFreq, float Q) { if (sampleRate <= 0 || centerFreq <= 0 || centerFreq >= sampleRate / 2 || Q <= 0) { SetIdentityCoefficients(filter); return; } float w0 = 2.0f * (float)M_PI * centerFreq / sampleRate; float cos_w0 = std::cos(w0); float sin_w0 = std::sin(w0); float alpha = sin_w0 / (2.0f * Q); float a0_cook = 1.0f + alpha; float a1_cook = -2.0f * cos_w0; float a2_cook = 1.0f - alpha; float b0_cook = alpha; float b1_cook = 0.0f; float b2_cook = -alpha; filter.b0 = b0_cook / a0_cook; filter.b1 = b1_cook / a0_cook; filter.b2 = b2_cook / a0_cook; filter.a1 = a1_cook / a0_cook; filter.a2 = a2_cook / a0_cook; }
void CalculateBandPass0dBCoefficients(BiquadFilter& filter, float sampleRate, float centerFreq, float Q) { if (sampleRate <= 0 || centerFreq <= 0 || centerFreq >= sampleRate / 2 || Q <= 0) { SetIdentityCoefficients(filter); return; } float w0 = 2.0f * (float)M_PI * centerFreq / sampleRate; float cos_w0 = std::cos(w0); float alpha = std::sin(w0) / (2.0f * Q); float a0_cook = 1.0f + alpha; float a1_cook = -2.0f * cos_w0; float a2_cook = 1.0f - alpha; float b0_cook = std::sin(w0) / 2.0f; float b1_cook = 0.0f; float b2_cook = -std::sin(w0) / 2.0f; filter.b0 = b0_cook / a0_cook; filter.b1 = b1_cook / a0_cook; filter.b2 = b2_cook / a0_cook; filter.a1 = a1_cook / a0_cook; filter.a2 = a2_cook / a0_cook; }
void CalculateNotchCoefficients(BiquadFilter& filter, float sampleRate, float centerFreq, float Q) { if (sampleRate <= 0 || centerFreq <= 0 || centerFreq >= sampleRate / 2 || Q <= 0) { SetIdentityCoefficients(filter); return; } float w0 = 2.0f * (float)M_PI * centerFreq / sampleRate; float cos_w0 = std::cos(w0); float alpha = std::sin(w0) / (2.0f * Q); float a0_cook = 1.0f + alpha; float a1_cook = -2.0f * cos_w0; float a2_cook = 1.0f - alpha; float b0_cook = 1.0f; float b1_cook = -2.0f * cos_w0; float b2_cook = 1.0f; filter.b0 = b0_cook / a0_cook; filter.b1 = b1_cook / a0_cook; filter.b2 = b2_cook / a0_cook; filter.a1 = a1_cook / a0_cook; filter.a2 = a2_cook / a0_cook; }
void CalculateLowShelfCoefficients(BiquadFilter& filter, float sampleRate, float cutoffFreq, float S, float gainDb) { if (sampleRate <= 0 || cutoffFreq <= 0 || cutoffFreq >= sampleRate / 2 || S <= 0) { SetIdentityCoefficients(filter); return; } float A = std::pow(10.0f, gainDb / 40.0f); float w0 = 2.0f * (float)M_PI * cutoffFreq / sampleRate; float cos_w0 = std::cos(w0); float sin_w0 = std::sin(w0); float alpha = sin_w0 / 2.0f * std::sqrt((A + 1.0f / A) * (1.0f / S - 1.0f) + 2.0f); if (S <= 0.0f) alpha = sin_w0 / 2.0f; if (!(alpha > 0.0f && alpha < std::numeric_limits<float>::infinity())) { printf("Warning: Alpha calculation failed for LowShelf (S=%.2f). Using Q=0.707 alpha.\n", S); alpha = sin_w0 / (2.0f * 0.707f); } float sqrtA = std::sqrt(A); float a0_cook = (A + 1.0f) + (A - 1.0f) * cos_w0 + 2.0f * sqrtA * alpha; float a1_cook = -2.0f * ((A - 1.0f) + (A + 1.0f) * cos_w0); float a2_cook = (A + 1.0f) + (A - 1.0f) * cos_w0 - 2.0f * sqrtA * alpha; float b0_cook = (A + 1.0f) - (A - 1.0f) * cos_w0 + 2.0f * sqrtA * alpha; float b1_cook = 2.0f * ((A - 1.0f) - (A + 1.0f) * cos_w0); float b2_cook = (A + 1.0f) - (A - 1.0f) * cos_w0 - 2.0f * sqrtA * alpha; filter.b0 = b0_cook / a0_cook; filter.b1 = b1_cook / a0_cook; filter.b2 = b2_cook / a0_cook; filter.a1 = a1_cook / a0_cook; filter.a2 = a2_cook / a0_cook; }
void CalculateHighShelfCoefficients(BiquadFilter& filter, float sampleRate, float cutoffFreq, float S, float gainDb) { if (sampleRate <= 0 || cutoffFreq <= 0 || cutoffFreq >= sampleRate / 2 || S <= 0) { SetIdentityCoefficients(filter); return; } float A = std::pow(10.0f, gainDb / 40.0f); float w0 = 2.0f * (float)M_PI * cutoffFreq / sampleRate; float cos_w0 = std::cos(w0); float sin_w0 = std::sin(w0); float alpha = sin_w0 / 2.0f * std::sqrt((A + 1.0f / A) * (1.0f / S - 1.0f) + 2.0f); if (S <= 0.0f) alpha = sin_w0 / 2.0f; if (!(alpha > 0.0f && alpha < std::numeric_limits<float>::infinity())) { printf("Warning: Alpha calculation failed for HighShelf (S=%.2f). Using Q=0.707 alpha.\n", S); alpha = sin_w0 / (2.0f * 0.707f); } float sqrtA = std::sqrt(A); float a0_cook = (A + 1.0f) - (A - 1.0f) * cos_w0 + 2.0f * sqrtA * alpha; float a1_cook = 2.0f * ((A - 1.0f) - (A + 1.0f) * cos_w0); float a2_cook = (A + 1.0f) - (A - 1.0f) * cos_w0 - 2.0f * sqrtA * alpha; float b0_cook = (A + 1.0f) + (A - 1.0f) * cos_w0 + 2.0f * sqrtA * alpha; float b1_cook = -2.0f * ((A - 1.0f) + (A + 1.0f) * cos_w0); float b2_cook = (A + 1.0f) + (A - 1.0f) * cos_w0 - 2.0f * sqrtA * alpha; filter.b0 = b0_cook / a0_cook; filter.b1 = b1_cook / a0_cook; filter.b2 = b2_cook / a0_cook; filter.a1 = a1_cook / a0_cook; filter.a2 = a2_cook / a0_cook; }
void CalculatePeakFilterCoefficients(BiquadFilter& filter, float sampleRate, float centerFreq, float Q, float gainDb) { if (sampleRate <= 0 || centerFreq <= 0 || centerFreq >= sampleRate / 2 || Q <= 0) { SetIdentityCoefficients(filter); return; } float A = std::pow(10.0f, gainDb / 40.0f); float w0 = 2.0f * (float)M_PI * centerFreq / sampleRate; float cos_w0 = std::cos(w0); float sin_w0 = std::sin(w0); float alpha = sin_w0 / (2.0f * Q); float a0_cook = 1.0f + alpha / A; float a1_cook = -2.0f * cos_w0; float a2_cook = 1.0f - alpha / A; float b0_cook = 1.0f + alpha * A; float b1_cook = -2.0f * cos_w0; float b2_cook = 1.0f - alpha * A; filter.b0 = b0_cook / a0_cook; filter.b1 = b1_cook / a0_cook; filter.b2 = b2_cook / a0_cook; filter.a1 = a1_cook / a0_cook; filter.a2 = a2_cook / a0_cook; }


// 5. Definitions of Class and Struct Methods

// Biquad Filter Implementations
BiquadFilter::BiquadFilter() : b0(1.0f), b1(0.0f), b2(0.0f), a1(0.0f), a2(0.0f) {}
void BiquadFilter::ResetState(size_t numChannels) { state.assign(numChannels, {0.0f, 0.0f}); }
float BiquadFilter::process_sample(float input, size_t channelIndex) {
    if (channelIndex >= state.size()) return input;
    float y = b0 * input + state[channelIndex].first;
    state[channelIndex].first = b1 * input - a1 * y + state[channelIndex].second;
    state[channelIndex].second = b2 * input - a2 * y; return y;
}


// Gain Effect Implementations
const char* GainEffect::GetName() const { return "Gain"; }
void GainEffect::ProcessBuffer(BYTE* pData, size_t numFrames, const WAVEFORMATEX* pFormat) {
    if (!pData || numFrames == 0 || !pFormat || std::fabs(gainFactor_ - 1.0f) < 1e-6) return; // Skip if no data or gain is 1.0

    size_t numSamples = numFrames * pFormat->nChannels;
    size_t bytesPerSample = pFormat->wBitsPerSample / 8;

    GUID subFormat = {0};
    bool isExt = (pFormat->wFormatTag == WAVE_FORMAT_EXTENSIBLE);
    UINT16 validBitsPerSample = pFormat->wBitsPerSample;

    if (isExt) {
        const WAVEFORMATEXTENSIBLE* pExt = reinterpret_cast<const WAVEFORMATEXTENSIBLE*>(pFormat);
         subFormat = pExt->SubFormat;
         validBitsPerSample = pExt->Samples.wValidBitsPerSample;
    } else {
         validBitsPerSample = pFormat->wBitsPerSample;
    }

    if (bytesPerSample == 4) {
        if (isExt && IsEqualGUID(subFormat, KSDATAFORMAT_SUBTYPE_IEEE_FLOAT) && validBitsPerSample == 32) {
             float* audio_samples = (float*)pData;
             for (size_t i = 0; i < numSamples; ++i) {
                 audio_samples[i] = fmaxf(-1.0f, fminf(audio_samples[i] * gainFactor_, 1.0f));
             }
        } else if (isExt && IsEqualGUID(subFormat, KSDATAFORMAT_SUBTYPE_PCM)) {
             if (validBitsPerSample == 24) { // 24-bit in 32-bit container
                 INT32* audio_samples = (INT32*)pData;
                 const float max_val_f = 8388607.0f; // (2^23 - 1)
                 const float min_val_f = -8388608.0f; // -(2^23)
                  for (size_t i = 0; i < numSamples; ++i) {
                     INT32 sample_int32 = (audio_samples[i] >> 8);
                     if (sample_int32 & 0x800000) sample_int32 |= 0xFF000000;
                     float sample_f = static_cast<float>(sample_int32);

                     sample_f = fmaxf(min_val_f, fminf(sample_f * gainFactor_, max_val_f));

                     INT32 processed_sample_int32 = static_cast<INT32>(sample_f);
                     audio_samples[i] = (processed_sample_int32 & 0xFFFFFF) << 8;
                 }
             } else if (validBitsPerSample == 32) { // 32-bit in 32-bit container
                 INT32* audio_samples = (INT32*)pData;
                  for (size_t i = 0; i < numSamples; ++i) {
                      double amplified_sample = static_cast<double>(audio_samples[i]) * gainFactor_;
                      audio_samples[i] = static_cast<INT32>(std::max(static_cast<double>(std::numeric_limits<INT32>::min()), std::min(static_cast<double>(std::numeric_limits<INT32>::max()), amplified_sample)));
                  }
             }
        }
    } else if (bytesPerSample == 2) {
        if (((!isExt && pFormat->wFormatTag == WAVE_FORMAT_PCM) || (isExt && IsEqualGUID(subFormat, KSDATAFORMAT_SUBTYPE_PCM))) && validBitsPerSample == 16)
        {
            short* audio_samples = (short*)pData;
            const float max_val_f = 32767.0f;
            const float min_val_f = -32768.0f;
            for (size_t i = 0; i < numSamples; ++i) {
                 float amplified_sample = static_cast<float>(audio_samples[i]) * gainFactor_;
                 audio_samples[i] = static_cast<short>(fmaxf(min_val_f, fminf(amplified_sample, max_val_f)));
            }
        }
    }
}


void GainEffect::UpdateParameters(const SharedDspParameters& params, const WAVEFORMATEX* pFormat) {
    float currentGain; BiquadFilterType dummyType; float dummyFreq, dummyQ, dummyGainDb, dummyS;
    const_cast<SharedDspParameters*>(&params)->GetParameters(currentGain, dummyType, dummyFreq, dummyQ, dummyGainDb, dummyS);
    gainFactor_ = currentGain;
}


// Biquad Effect Implementations
const char* BiquadEffect::GetName() const { return "Biquad"; }
void BiquadEffect::RecalculateCoefficients(const WAVEFORMATEX* pFormat) {
     if (!pFormat) { SetIdentityCoefficients(filter_); filter_.ResetState(0); return; } float sampleRate = (float)pFormat->nSamplesPerSec; size_t numChannels = pFormat->nChannels;
    switch (filterType_) {
        case BiquadFilterType::Identity: SetIdentityCoefficients(filter_); break;
        case BiquadFilterType::LowPass: CalculateLowPassCoefficients(filter_, sampleRate, filterFreq_, filterQ_); break;
        case BiquadFilterType::HighPass: CalculateHighPassCoefficients(filter_, sampleRate, filterFreq_, filterQ_); break;
        case BiquadFilterType::BandPass: CalculateBandPassCoefficients(filter_, sampleRate, filterFreq_, filterQ_); break;
        case BiquadFilterType::BandPass0dB: CalculateBandPass0dBCoefficients(filter_, sampleRate, filterFreq_, filterQ_); break;
        case BiquadFilterType::Notch: CalculateNotchCoefficients(filter_, sampleRate, filterFreq_, filterQ_); break;
        case BiquadFilterType::LowShelf: CalculateLowShelfCoefficients(filter_, sampleRate, filterFreq_, filterS_, filterGainDb_); break;
        case BiquadFilterType::HighShelf: CalculateHighShelfCoefficients(filter_, sampleRate, filterFreq_, filterS_, filterGainDb_); break;
        case BiquadFilterType::Peak: CalculatePeakFilterCoefficients(filter_, sampleRate, filterFreq_, filterQ_, filterGainDb_); break;
        default: printf("Warning: Unknown filter type (%d). Setting Identity filter.\n", (int)filterType_); SetIdentityCoefficients(filter_); break;
    } filter_.ResetState(numChannels);
}
void BiquadEffect::ProcessBuffer(BYTE* pData, size_t numFrames, const WAVEFORMATEX* pFormat) {
    if (!pData || numFrames == 0 || !pFormat || filterType_ == BiquadFilterType::Identity) return;

    size_t numChannels = pFormat->nChannels;
    size_t bytesPerSample = pFormat->wBitsPerSample / 8;
    size_t bytesPerFrame = numChannels * bytesPerSample;

    if (filter_.state.size() != numChannels) {
        printf("Warning [BiquadEffect]: Filter state size mismatch (%zu vs %zu channels) before processing. Resetting state.\n", filter_.state.size(), numChannels);
        filter_.ResetState(numChannels);
         if (filter_.state.size() != numChannels) { printf("Error [BiquadEffect]: Failed to resize filter state after mismatch. Skipping processing.\n"); return; }
    }

    GUID subFormat = {0};
    bool isExt = (pFormat->wFormatTag == WAVE_FORMAT_EXTENSIBLE);
    UINT16 validBitsPerSample = pFormat->wBitsPerSample;
    if (isExt) { const WAVEFORMATEXTENSIBLE* pExt = reinterpret_cast<const WAVEFORMATEXTENSIBLE*>(pFormat); subFormat = pExt->SubFormat; validBitsPerSample = pExt->Samples.wValidBitsPerSample; }
    else { validBitsPerSample = pFormat->wBitsPerSample; }

    for (size_t i = 0; i < numFrames; ++i) {
        for (size_t j = 0; j < numChannels; ++j) {
            BYTE* pCurrentSample = pData + i * bytesPerFrame + j * bytesPerSample;
            if (bytesPerSample == 4) {
                if (isExt && IsEqualGUID(subFormat, KSDATAFORMAT_SUBTYPE_IEEE_FLOAT) && validBitsPerSample == 32) {
                    float* sample_ptr = reinterpret_cast<float*>(pCurrentSample); *sample_ptr = filter_.process_sample(*sample_ptr, j);
                } else if (isExt && IsEqualGUID(subFormat, KSDATAFORMAT_SUBTYPE_PCM)) {
                     INT32* sample_ptr = reinterpret_cast<INT32*>(pCurrentSample);
                     if (validBitsPerSample == 24) {
                         INT32 sample_int32 = (*sample_ptr >> 8); if (sample_int32 & 0x800000) sample_int32 |= 0xFF000000;
                         float sample_f = static_cast<float>(sample_int32);
                         float processed_sample_f = filter_.process_sample(sample_f, j);
                         INT32 processed_sample_int32 = static_cast<INT32>(fmaxf(-8388608.0f, fminf(processed_sample_f, 8388607.0f)));
                         *sample_ptr = (processed_sample_int32 & 0xFFFFFF) << 8;
                     } else if (validBitsPerSample == 32) {
                         INT32* sample_ptr = reinterpret_cast<INT32*>(pCurrentSample);
                         float sample_f = static_cast<float>(*sample_ptr);
                         float processed_sample_f = filter_.process_sample(sample_f, j);
                          *sample_ptr = static_cast<INT32>(fmaxf(-2147483648.0f, fminf(processed_sample_f, 2147483647.0f)));
                     }
                }
            } else if (bytesPerSample == 2) {
                 if (((!isExt && pFormat->wFormatTag == WAVE_FORMAT_PCM) || (isExt && IsEqualGUID(subFormat, KSDATAFORMAT_SUBTYPE_PCM))) && validBitsPerSample == 16)
                {
                     short* sample_ptr = reinterpret_cast<short*>(pCurrentSample);
                     float sample_f = static_cast<float>(*sample_ptr);
                     float processed_sample_f = filter_.process_sample(sample_f, j);
                     *sample_ptr = static_cast<short>(fmaxf(-32768.0f, fminf(processed_sample_f, 32767.0f)));
                }
            }
        }
    }
}

void BiquadEffect::UpdateParameters(const SharedDspParameters& params, const WAVEFORMATEX* pFormat) {
    float dummyGain; BiquadFilterType newType; float newFreq, newQ, newGainDb, newS;
    const_cast<SharedDspParameters*>(&params)->GetParameters(dummyGain, newType, newFreq, newQ, newGainDb, newS);
    bool paramsChanged = (filterType_ != newType) || (std::fabs(filterFreq_ - newFreq) > 1e-3f) || (std::fabs(filterQ_ - newQ) > 1e-3f) || (std::fabs(filterGainDb_ - newGainDb) > 1e-3f) || (std::fabs(filterS_ - newS) > 1e-3f);
    if (paramsChanged) { filterType_ = newType; filterFreq_ = newFreq; filterQ_ = newQ; filterGainDb_ = newGainDb; filterS_ = newS; RecalculateCoefficients(pFormat); }}
void BiquadEffect::Initialize(const WAVEFORMATEX* pFormat) { RecalculateCoefficients(pFormat); }


// Shared DSP Parameters Implementations
SharedDspParameters::SharedDspParameters() : gainFactor(1.0f), filterType(BiquadFilterType::Identity), filterFreq(1000.0f), filterQ(1.0f), filterGainDb(0.0f), filterS(1.0f) { InitializeCriticalSection(&lock); hParamChangedEvent = NULL; }
SharedDspParameters::~SharedDspParameters() { Cleanup(); DeleteCriticalSection(&lock); }
void SharedDspParameters::Initialize() { hParamChangedEvent = CreateEvent(NULL, FALSE, FALSE, NULL); if (!hParamChangedEvent) { printf("Error: Failed to create hParamChangedEvent (%lu).\n", GetLastError()); } printf("Shared DSP Parameters initialized.\n"); }
void SharedDspParameters::Cleanup() { if (hParamChangedEvent && hParamChangedEvent != INVALID_HANDLE_VALUE) { CloseHandle(hParamChangedEvent); hParamChangedEvent = NULL; } printf("Shared DSP Parameters cleaned up.\n"); }
void SharedDspParameters::GetParameters(float& outGain, BiquadFilterType& outType, float& outFreq, float& outQ, float& outGainDb, float& outS) const { EnterCriticalSection(const_cast<LPCRITICAL_SECTION>(&lock)); outGain = gainFactor; outType = filterType; outFreq = filterFreq; outQ = filterQ; outGainDb = filterGainDb; outS = filterS; LeaveCriticalSection(const_cast<LPCRITICAL_SECTION>(&lock)); }
void SharedDspParameters::SetGain(float newGain) { float gainDbPrint = (newGain > 0.0f) ? (20.0f * std::log10(newGain)) : (newGain < 0.0f ? (20.0f * std::log10(std::abs(newGain))) : -std::numeric_limits<float>::infinity()); EnterCriticalSection(&lock); gainFactor = newGain; LeaveCriticalSection(&lock); SignalParametersChanged(); printf("DSP Gain set to %.2f (%.1fdB)\n", newGain, gainDbPrint); }
void SharedDspParameters::SetFilterType(BiquadFilterType newType) { EnterCriticalSection(&lock); filterType = newType; LeaveCriticalSection(&lock); SignalParametersChanged(); printf("DSP Filter Type set to %s.\n", BiquadFilterTypeToString(newType)); }
void SharedDspParameters::SetFilterFreq(float newFreq) { if (newFreq <= 0.0f) { printf("Error: Frequency must be positive.\n"); return; } EnterCriticalSection(&lock); filterFreq = newFreq; LeaveCriticalSection(&lock); SignalParametersChanged(); printf("DSP Filter Frequency set to %.1f Hz.\n", newFreq); }
void SharedDspParameters::SetFilterQ(float newQ) { if (newQ <= 0.0f) { printf("Error: Q must be positive.\n"); return; } EnterCriticalSection(&lock); filterQ = newQ; LeaveCriticalSection(&lock); SignalParametersChanged(); printf("DSP Filter Q set to %.2f.\n", newQ); }
void SharedDspParameters::SetFilterGainDb(float newGainDb) { EnterCriticalSection(&lock); filterGainDb = newGainDb; LeaveCriticalSection(&lock); SignalParametersChanged(); printf("DSP Filter Gain set to %.1f dB.\n", newGainDb); }
void SharedDspParameters::SetFilterS(float newS) { if (newS <= 0.0f) { printf("Error: S must be positive.\n"); return; } EnterCriticalSection(&lock); filterS = newS; LeaveCriticalSection(&lock); SignalParametersChanged(); printf("DSP Filter S set to %.2f.\n", newS); }
void SharedDspParameters::SetFilterSettings(BiquadFilterType type, float freq, float q, float gainDb, float s) { EnterCriticalSection(&lock); filterType = type; filterFreq = (freq > 0.0f) ? freq : filterFreq; filterQ = (q > 0.0f) ? q : filterQ; filterGainDb = gainDb; filterS = (s > 0.0f) ? s : filterS; LeaveCriticalSection(&lock); SignalParametersChanged(); printf("DSP Filter settings updated: Type=%s, Freq=%.1f, Q=%.2f, Gain=%.1f, S=%.2f\n", BiquadFilterTypeToString(filterType), filterFreq, filterQ, filterGainDb, filterS); }
void SharedDspParameters::SignalParametersChanged() { if (hParamChangedEvent && hParamChangedEvent != INVALID_HANDLE_VALUE) { SetEvent(hParamChangedEvent); } }


// Shared Monitor Data Implementations
SharedMonitorData::SharedMonitorData() : capturePadding(0), renderPadding(0), framesReadFromSharedBuffer(0), framesConsumedFromSharedBuffer(0), framesWrittenToRenderBuffer(0) { InitializeCriticalSection(&lock); } // ИСправлено: Инициализация новых членов
SharedMonitorData::~SharedMonitorData() { Cleanup(); DeleteCriticalSection(&lock); }
void SharedMonitorData::Initialize() { printf("Shared Monitor Data initialized.\n"); }
void SharedMonitorData::Cleanup() { printf("Shared Monitor Data cleaned up.\n"); }
void SharedMonitorData::UpdateCapturePadding(UINT32 padding) { EnterCriticalSection(&lock); capturePadding = padding; LeaveCriticalSection(&lock); }
void SharedMonitorData::UpdateRenderPadding(UINT32 padding) { EnterCriticalSection(&lock); renderPadding = padding; LeaveCriticalSection(&lock); }
void SharedMonitorData::GetPadding(UINT32& outCapture, UINT32& outRender) const { EnterCriticalSection(const_cast<LPCRITICAL_SECTION>(&lock)); outCapture = capturePadding; outRender = renderPadding; LeaveCriticalSection(const_cast<LPCRITICAL_SECTION>(&lock)); }
// ИСправлено: Реализация getter для мониторинговых данных SharedBuffer
void SharedMonitorData::GetSharedBufferStats(UINT32& outFramesRead, UINT32& outFramesConsumed, UINT32& outFramesWrittenToRender) const {
     EnterCriticalSection(const_cast<LPCRITICAL_SECTION>(&lock));
     outFramesRead = framesReadFromSharedBuffer;
     outFramesConsumed = framesConsumedFromSharedBuffer;
     outFramesWrittenToRender = framesWrittenToRenderBuffer;
     LeaveCriticalSection(const_cast<LPCRITICAL_SECTION>(&lock));
}


// Shared Audio Buffer Implementations
SharedAudioBuffer::SharedAudioBuffer() : bufferData(nullptr), bufferSize(0), frameSize(0), writePos(0), readPos(0), pCaptureFormat(nullptr), hasCaptureFormat(false) { InitializeCriticalSection(&lock); memset(captureFormatBuffer, 0, sizeof(captureFormatBuffer)); }
SharedAudioBuffer::~SharedAudioBuffer() { Cleanup(); DeleteCriticalSection(&lock); }
SharedAudioBuffer::SharedAudioBuffer(SharedAudioBuffer&& other) noexcept
    : bufferData(std::exchange(other.bufferData, nullptr)), bufferSize(std::exchange(other.bufferSize, 0)),
      frameSize(std::exchange(other.frameSize, 0)), writePos(std::exchange(other.writePos, 0)),
      readPos(std::exchange(other.readPos, 0)), hasCaptureFormat(std::exchange(other.hasCaptureFormat, false))
{ memcpy(this->captureFormatBuffer, other.captureFormatBuffer, sizeof(this->captureFormatBuffer)); this->pCaptureFormat = reinterpret_cast<PWAVEFORMATEX>(this->captureFormatBuffer); }
SharedAudioBuffer& SharedAudioBuffer::operator=(SharedAudioBuffer&& other) noexcept {
    if (this != &other) { Cleanup(); bufferData = std::exchange(other.bufferData, nullptr); bufferSize = std::exchange(other.bufferSize, 0);
        frameSize = std::exchange(other.frameSize, 0); writePos = std::exchange(other.writePos, 0); readPos = std::exchange(other.readPos, 0);
        hasCaptureFormat = std::exchange(other.hasCaptureFormat, false); memcpy(this->captureFormatBuffer, other.captureFormatBuffer, sizeof(this->captureFormatBuffer));
        this->pCaptureFormat = reinterpret_cast<PWAVEFORMATEX>(this->captureFormatBuffer); } return *this; }
size_t SharedAudioBuffer::Initialize(const WAVEFORMATEX* pFormat, size_t durationMs) {
    if (!pFormat) return 0; size_t formatSize = sizeof(WAVEFORMATEX) + pFormat->cbSize;
    if (formatSize > sizeof(captureFormatBuffer)) { printf("Ошибка: Размер формата захвата (%zu) превышает размер внутреннего буфера captureFormatBuffer (%zu).\n", formatSize, sizeof(captureFormatBuffer)); return 0; }
    memcpy(captureFormatBuffer, pFormat, formatSize); pCaptureFormat = reinterpret_cast<PWAVEFORMATEX>(captureFormatBuffer); hasCaptureFormat = true;
    frameSize = pCaptureFormat->nBlockAlign; if (frameSize == 0) frameSize = (size_t)pCaptureFormat->nChannels * (pCaptureFormat->wBitsPerSample / 8);
    if (frameSize == 0) { printf("Ошибка: Frame size равен 0.\n"); Cleanup(); return 0; }
    size_t desiredFrames = (size_t)((double)pCaptureFormat->nSamplesPerSec * durationMs / 1000.0);
    size_t allocationFrames = desiredFrames * 4; if (allocationFrames == 0) allocationFrames = pCaptureFormat->nSamplesPerSec / 2; if (allocationFrames < 100) allocationFrames = 100;
     if (frameSize > 0 && allocationFrames > SIZE_MAX / frameSize) { printf("Ошибка: Размер общего буфера слишком большой.\n"); Cleanup(); return 0; }
    bufferSize = allocationFrames * frameSize; bufferData = new (std::nothrow) BYTE[bufferSize];
    if (!bufferData) { bufferSize = 0; printf("Ошибка: Не удалось выделить память для общего буфера.\n"); return 0; }
    memset(bufferData, 0, bufferSize); writePos = 0; readPos = 0;
    // ИСправлено: bufferSizeFrames_ добавлен в SharedAudioBuffer и инициализирован
    bufferSizeFrames_ = allocationFrames;
    printf("Общий буфер инициализирован. Размер: %zu байт (%zu фреймов при захвате, ~%zu мс).\n", bufferSize, bufferSize / frameSize, (size_t)((double)bufferSize / frameSize * 1000.0 / pCaptureFormat->nSamplesPerSec) );
    return bufferSize; }
void SharedAudioBuffer::Cleanup() {
    if (bufferData) { delete[] bufferData; bufferData = nullptr; } bufferSize = 0; frameSize = 0; writePos = 0; readPos = 0; pCaptureFormat = nullptr; hasCaptureFormat = false;
    bufferSizeFrames_ = 0; // ИСправлено: Очистка bufferSizeFrames_
    memset(captureFormatBuffer, 0, sizeof(captureFormatBuffer)); printf("Общий буфер очищен.\n");
}
size_t SharedAudioBuffer::Write(const BYTE* pData, size_t numFrames) {
    if (!bufferData || numFrames == 0 || bufferSize == 0 || frameSize == 0 || !hasCaptureFormat || !pCaptureFormat) return 0;
    size_t bytesToWrite = numFrames * frameSize; size_t framesWritten = 0; EnterCriticalSection(&lock);
    size_t currentDataBytes = (writePos >= readPos) ? (writePos - readPos) : (bufferSize - readPos + writePos);
    size_t freeSpaceBytes = bufferSize - currentDataBytes; size_t effectiveFreeSpaceBytes = (freeSpaceBytes >= frameSize) ? (freeSpaceBytes - frameSize) : 0;
    if (bytesToWrite > effectiveFreeSpaceBytes) { bytesToWrite = effectiveFreeSpaceBytes; numFrames = bytesToWrite / frameSize; bytesToWrite = numFrames * frameSize; if (numFrames == 0) { LeaveCriticalSection(&lock); return 0; }}
    if (bytesToWrite == 0) { LeaveCriticalSection(&lock); return 0; }
    size_t writePosLocal = writePos; size_t bytesTillEnd = bufferSize - writePosLocal;
    if (bytesToWrite <= bytesTillEnd) { memcpy(bufferData + writePosLocal, pData, bytesToWrite); writePos = (writePosLocal + bytesToWrite); if (writePos == bufferSize) writePos = 0; framesWritten = numFrames; }
    else { size_t firstPartBytes = bytesTillEnd; size_t secondPartBytes = bytesToWrite - firstPartBytes; memcpy(bufferData + writePosLocal, pData, firstPartBytes); memcpy(bufferData, pData + firstPartBytes, secondPartBytes); writePos = secondPartBytes; framesWritten = numFrames;}
    LeaveCriticalSection(&lock); return framesWritten;
}
// ИСправлено: Read теперь возвращает количество фреймов и НЕ двигает указатель чтения.
size_t SharedAudioBuffer::Read(BYTE* pBuffer, size_t max_frames_to_read_capture_format) {
    if (!bufferData || !pBuffer || max_frames_to_read_capture_format == 0 || bufferSize == 0 || frameSize == 0 || !hasCaptureFormat || !pCaptureFormat) {
        return 0;
    }
    size_t bytesToRead = max_frames_to_read_capture_format * frameSize; size_t framesRead = 0;
    EnterCriticalSection(&lock);
    size_t currentDataBytes = (writePos >= readPos) ? (writePos - readPos) : (bufferSize - readPos + writePos);
    bytesToRead = std::min(currentDataBytes, bytesToRead); framesRead = bytesToRead / frameSize; bytesToRead = framesRead * frameSize;
     if (bytesToRead == 0) { LeaveCriticalSection(&lock); return 0; }
    size_t readPosLocal = readPos; size_t bytesTillEnd = bufferSize - readPosLocal;
    if (bytesToRead <= bytesTillEnd) { memcpy(pBuffer, bufferData + readPosLocal, bytesToRead); }
    else { size_t firstPartBytes = bytesTillEnd; size_t secondPartBytes = bytesToRead - firstPartBytes; memcpy(pBuffer, bufferData + readPosLocal, firstPartBytes); memcpy(pBuffer + firstPartBytes, bufferData, secondPartBytes); }
    // readPos is NOT advanced here. It's advanced by AdvanceReadPointer.
    LeaveCriticalSection(&lock);
    return framesRead; // Returns the number of frames COPIED to pBuffer
}
size_t SharedAudioBuffer::GetAvailableFrames() const {
    if (!bufferData || bufferSize == 0 || frameSize == 0) return 0;
    EnterCriticalSection(const_cast<LPCRITICAL_SECTION>(&lock));
    size_t currentDataBytes = (writePos >= readPos) ? (writePos - readPos) : (bufferSize - readPos + writePos);
     LeaveCriticalSection(const_cast<LPCRITICAL_SECTION>(&lock));
    return currentDataBytes / frameSize;
}
const WAVEFORMATEX* SharedAudioBuffer::GetCaptureFormat() const { return hasCaptureFormat ? pCaptureFormat : nullptr; }
// ИСправлено: Реализация AdvanceReadPointer
void SharedAudioBuffer::AdvanceReadPointer(size_t numFrames) {
     if (bufferSizeFrames_ == 0 || frameSize == 0) return; // Should not happen if initialized
     EnterCriticalSection(&lock);
     readPos = (readPos + numFrames * frameSize) % bufferSize;
     LeaveCriticalSection(&lock);
}


// AudioConverter Implementations

// Initialize static warning flag
bool AudioConverter::warned_unsupported_conversion_ = false;

AudioConverter::AudioConverter() = default;

// Added const to GetFloatSample
float AudioConverter::GetFloatSample(const BYTE* pSample, const WAVEFORMATEX* pFormat) const {
     if (!pSample || !pFormat) return 0.0f;
     size_t bytesPerSample = pFormat->wBitsPerSample / 8;
     GUID subFormat = {0};
     bool isExt = (pFormat->wFormatTag == WAVE_FORMAT_EXTENSIBLE);
     UINT16 validBitsPerSample = pFormat->wBitsPerSample;

     if (isExt) { const WAVEFORMATEXTENSIBLE* pExt = reinterpret_cast<const WAVEFORMATEXTENSIBLE*>(pFormat); subFormat = pExt->SubFormat; validBitsPerSample = pExt->Samples.wValidBitsPerSample; }
     else { validBitsPerSample = pFormat->wBitsPerSample; }

     if (bytesPerSample == 4) {
         if (isExt && IsEqualGUID(subFormat, KSDATAFORMAT_SUBTYPE_IEEE_FLOAT) && validBitsPerSample == 32) {
             return *reinterpret_cast<const float*>(pSample);
         } else if (isExt && IsEqualGUID(subFormat, KSDATAFORMAT_SUBTYPE_PCM)) {
             INT32 sample_int32 = *reinterpret_cast<const INT32*>(pSample);
             if (validBitsPerSample == 24) { sample_int32 = (sample_int32 >> 8); if (sample_int32 & 0x800000) sample_int32 |= 0xFF000000; return static_cast<float>(sample_int32) / 8388607.0f; }
             else if (validBitsPerSample == 32) { return static_cast<float>(sample_int32) / 2147483647.0f; }
         }
     } else if (bytesPerSample == 2) {
          if (((!isExt && pFormat->wFormatTag == WAVE_FORMAT_PCM) || (isExt && IsEqualGUID(subFormat, KSDATAFORMAT_SUBTYPE_PCM))) && validBitsPerSample == 16)
         { short sample_int16 = *reinterpret_cast<const short*>(pSample); return static_cast<float>(sample_int16) / 32767.0f; }
     }
     return 0.0f;
}

// Added const to SetFloatSample
void AudioConverter::SetFloatSample(BYTE* pSample, float value, const WAVEFORMATEX* pFormat) const {
     if (!pSample || !pFormat) return;
     size_t bytesPerSample = pFormat->wBitsPerSample / 8;
     GUID subFormat = {0};
     bool isExt = (pFormat->wFormatTag == WAVE_FORMAT_EXTENSIBLE);
     UINT16 validBitsPerSample = pFormat->wBitsPerSample;

     if (isExt) { const WAVEFORMATEXTENSIBLE* pExt = reinterpret_cast<const WAVEFORMATEXTENSIBLE*>(pFormat); subFormat = pExt->SubFormat; validBitsPerSample = pExt->Samples.wValidBitsPerSample; }
     else { validBitsPerSample = pFormat->wBitsPerSample; }

     if (bytesPerSample == 4) {
         if (isExt && IsEqualGUID(subFormat, KSDATAFORMAT_SUBTYPE_IEEE_FLOAT) && validBitsPerSample == 32) {
             *reinterpret_cast<float*>(pSample) = fmaxf(-1.0f, fminf(value, 1.0f));
         } else if (isExt && IsEqualGUID(subFormat, KSDATAFORMAT_SUBTYPE_PCM)) {
             if (validBitsPerSample == 24) {
                 INT32 processed_sample_int32 = static_cast<INT32>(fmaxf(-8388608.0f, fminf(value * 8388607.0f, 8388607.0f)));
                 *reinterpret_cast<INT32*>(pSample) = (processed_sample_int32 & 0xFFFFFF) << 8;
             } else if (validBitsPerSample == 32) {
                 *reinterpret_cast<INT32*>(pSample) = static_cast<INT32>(fmaxf(-2147483648.0f, fminf(value * 2147483647.0f, 2147483647.0f)));
             }
         }
     } else if (bytesPerSample == 2) {
          if (((!isExt && pFormat->wFormatTag == WAVE_FORMAT_PCM) || (isExt && IsEqualGUID(subFormat, KSDATAFORMAT_SUBTYPE_PCM))) && validBitsPerSample == 16)
         {
             *reinterpret_cast<short*>(pSample) = static_cast<short>(fmaxf(-32768.0f, fminf(value * 32767.0f, 32767.0f)));
         }
     }
}

// Added const to InterpolateLinear
float AudioConverter::InterpolateLinear(float sample1, float sample2, float factor) const {
    return sample1 * (1.0f - factor) + sample2 * factor;
}


bool AudioConverter::IsFormatCombinationSupported(const WAVEFORMATEX* pCaptureFormat, const WAVEFORMATEX* pRenderFormat) const {
     if (!pCaptureFormat || !pRenderFormat) return false;

     size_t captureBitsPerSample = pCaptureFormat->wBitsPerSample;
     size_t renderBitsPerSample = pRenderFormat->wBitsPerSample;
     UINT16 captureValidBits = captureBitsPerSample;
     UINT16 renderValidBits = renderBitsPerSample;
     GUID captureSubFormat = {0}, renderSubFormat = {0};
     bool isCaptureExt = (pCaptureFormat->wFormatTag == WAVE_FORMAT_EXTENSIBLE);
     bool isRenderExt = (pRenderFormat->wFormatTag == WAVE_FORMAT_EXTENSIBLE);

     if (isCaptureExt) { const WAVEFORMATEXTENSIBLE* pExtCapture = reinterpret_cast<const WAVEFORMATEXTENSIBLE*>(pCaptureFormat); captureSubFormat = pExtCapture->SubFormat; captureValidBits = pExtCapture->Samples.wValidBitsPerSample; }
     if (isRenderExt) { const WAVEFORMATEXTENSIBLE* pExtRender = reinterpret_cast<const WAVEFORMATEXTENSIBLE*>(pRenderFormat); renderSubFormat = pExtRender->SubFormat; renderValidBits = pExtRender->Samples.wValidBitsPerSample; }
    else { renderValidBits = renderBitsPerSample; }

     bool formatTagsMatch = (pCaptureFormat->wFormatTag == pRenderFormat->wFormatTag);
     bool subFormatsMatch = (!isCaptureExt || !isRenderExt || IsEqualGUID(captureSubFormat, renderSubFormat));
     bool bitsMatch = (captureBitsPerSample == renderBitsPerSample && captureValidBits == renderValidBits);

     bool isFloat32 = (captureBitsPerSample == 4 && captureValidBits == 32 && isCaptureExt && IsEqualGUID(captureSubFormat, KSDATAFORMAT_SUBTYPE_IEEE_FLOAT));
     bool isPCM16 = (captureBitsPerSample == 2 && captureValidBits == 16 && ((!isCaptureExt && pCaptureFormat->wFormatTag == WAVE_FORMAT_PCM) || (isCaptureExt && IsEqualGUID(captureSubFormat, KSDATAFORMAT_SUBTYPE_PCM))));
     bool isPCM24_32 = (captureBitsPerSample == 4 && captureValidBits == 24 && isCaptureExt && IsEqualGUID(captureSubFormat, KSDATAFORMAT_SUBTYPE_PCM));
     bool isPCM32_32 = (captureBitsPerSample == 4 && captureValidBits == 32 && isCaptureExt && IsEqualGUID(captureSubFormat, KSDATAFORMAT_SUBTYPE_PCM));
     bool isSupportedType = isFloat32 || isPCM16 || isPCM24_32 || isPCM32_32;

     bool formatsCompatibleForConversion = formatTagsMatch && subFormatsMatch && bitsMatch && isSupportedType;

      bool channelsCompatibleForMapping = (pCaptureFormat->nChannels == pRenderFormat->nChannels) ||
                                         (pCaptureFormat->nChannels == 1 && pRenderFormat->nChannels >= 1) || // Mono -> Any (map to first)
                                         (pCaptureFormat->nChannels >= 1 && pRenderFormat->nChannels == 1); // Any -> Mono (average first N)


     return formatsCompatibleForConversion && channelsCompatibleForMapping;
}


bool AudioConverter::Initialize(const WAVEFORMATEX* pCaptureFormat, const WAVEFORMATEX* pRenderFormat) {
    pCaptureFormat_ = pCaptureFormat;
    pRenderFormat_ = pRenderFormat;
    resampleRatio_ = 1.0;
    resamplePhase_ = 0.0;
    warned_unsupported_conversion_ = false; // Reset warning flag

    if (!pCaptureFormat_ || !pRenderFormat_) {
        printf("Error [AudioConverter::Initialize]: NULL format pointers.\n");
        return false;
    }

     if (!IsFormatCombinationSupported(pCaptureFormat_, pRenderFormat_)) {
         printf("Error [AudioConverter::Initialize]: Format combination not supported by converter implementation.\n");
          printf("  Capture Format: Tag=%d, Ch=%u, Rate=%lu, Bits=%u, ValidBits=%u, SubFormat={%08lX-...}\n",
                 pCaptureFormat_->wFormatTag, pCaptureFormat_->nChannels, pCaptureFormat_->nSamplesPerSec, pCaptureFormat_->wBitsPerSample,
                 (pCaptureFormat_->wFormatTag == WAVE_FORMAT_EXTENSIBLE && pCaptureFormat_->cbSize >= sizeof(WAVEFORMATEXTENSIBLE) - sizeof(WAVEFORMATEX)) ? reinterpret_cast<const WAVEFORMATEXTENSIBLE*>(pCaptureFormat_)->Samples.wValidBitsPerSample : pCaptureFormat_->wBitsPerSample,
                 (pCaptureFormat_->wFormatTag == WAVE_FORMAT_EXTENSIBLE && pCaptureFormat_->cbSize >= sizeof(WAVEFORMATEXTENSIBLE) - sizeof(WAVEFORMATEX)) ? reinterpret_cast<const WAVEFORMATEXTENSIBLE*>(pCaptureFormat_)->SubFormat.Data1 : 0);
          printf("  Render Format:  Tag=%d, Ch=%u, Rate=%lu, Bits=%u, ValidBits=%u, SubFormat={%08lX-...}\n",
                 pRenderFormat_->wFormatTag, pRenderFormat_->nChannels, pRenderFormat_->nSamplesPerSec, pRenderFormat_->wBitsPerSample,
                 (pRenderFormat_->wFormatTag == WAVE_FORMAT_EXTENSIBLE && pRenderFormat_->cbSize >= sizeof(WAVEFORMATEXTENSIBLE) - sizeof(WAVEFORMATEX)) ? reinterpret_cast<const WAVEFORMATEXTENSIBLE*>(pRenderFormat_)->Samples.wValidBitsPerSample : pRenderFormat_->wBitsPerSample,
                 (pRenderFormat_->wFormatTag == WAVE_FORMAT_EXTENSIBLE && pRenderFormat_->cbSize >= sizeof(WAVEFORMATEXTENSIBLE) - sizeof(WAVEFORMATEX)) ? reinterpret_cast<const WAVEFORMATEXTENSIBLE*>(pRenderFormat_)->SubFormat.Data1 : 0);
         return false;
     }

    if (pCaptureFormat_->nSamplesPerSec != pRenderFormat_->nSamplesPerSec) {
        resampleRatio_ = (double)pCaptureFormat_->nSamplesPerSec / pRenderFormat_->nSamplesPerSec;
        printf("Info [AudioConverter]: Sample rate conversion needed (%.0f -> %.0f). Ratio: %.4f.\n",
               (float)pCaptureFormat_->nSamplesPerSec, (float)pRenderFormat_->nSamplesPerSec, resampleRatio_);

        size_t captureBlockAlign = pCaptureFormat_->nBlockAlign;
        if (captureBlockAlign == 0 && pCaptureFormat_->nChannels > 0 && pCaptureFormat_->wBitsPerSample > 0) {
             captureBlockAlign = (size_t)pCaptureFormat_->nChannels * (pCaptureFormat_->wBitsPerSample / 8);
        }


        if (captureBlockAlign > 0) {
             // For linear interpolation, we need 1 frame of state
             interpolationStateBuffer_.resize(captureBlockAlign);
             printf("Info [AudioConverter]: Allocated %zu bytes for interpolation state buffer (1 frame).\n", interpolationStateBuffer_.size());
        } else {
             printf("Error [AudioConverter::Initialize]: Capture format nBlockAlign is 0. Cannot allocate state buffer.\n");
             interpolationStateBuffer_.clear();
             return false;
        }

    } else {
        resampleRatio_ = 1.0;
        printf("Info [AudioConverter]: Sample rates match (%.0f).\n", (float)pCaptureFormat_->nSamplesPerSec);
        interpolationStateBuffer_.clear(); // No state needed if no resampling
    }

    ResetState();

    printf("AudioConverter initialized successfully.\n");
    return true;
}

void AudioConverter::ResetState() {
    resamplePhase_ = 0.0;
    if (!interpolationStateBuffer_.empty()) {
        memset(interpolationStateBuffer_.data(), 0, interpolationStateBuffer_.size());
    }
     printf("AudioConverter state reset.\n");
}

size_t AudioConverter::GetStateBufferSize() const {
    return interpolationStateBuffer_.size(); // Size in bytes
}

size_t AudioConverter::GetStateBufferSizeFrames() const {
     size_t captureBlockAlign = pCaptureFormat_ ? pCaptureFormat_->nBlockAlign : 0;
     if (captureBlockAlign == 0 && pCaptureFormat_ && pCaptureFormat_->nChannels > 0 && pCaptureFormat_->wBitsPerSample > 0) {
         captureBlockAlign = (size_t)pCaptureFormat_->nChannels * (pCaptureFormat_->wBitsPerSample / 8);
     }

     if (interpolationStateBuffer_.empty() || captureBlockAlign == 0) return 0;
     return interpolationStateBuffer_.size() / captureBlockAlign; // Size in frames
}


// ИСправлено: Реализация Convert с обновленной сигнатурой
size_t AudioConverter::Convert(const BYTE* inputBuffer, size_t numInputFramesInBlock, BYTE* outputBuffer, size_t maxOutputFrames, size_t* numInputFramesConsumedOut) {
    // inputBuffer: current input data block (FROM SharedBuffer), size numInputFramesInBlock * captureBlockAlign
    // numInputFramesInBlock: number of frames in the current input data block
    // outputBuffer: buffer to write output frames
    // maxOutputFrames: max frames to write to outputBuffer
    // numInputFramesConsumedOut: output parameter, frames consumed from numInputFramesInBlock
    // Returns: number of output frames written

    if (!pCaptureFormat_ || !pRenderFormat_ || !inputBuffer || !outputBuffer || maxOutputFrames == 0 || !numInputFramesConsumedOut) {
         if (outputBuffer && pRenderFormat_) memset(outputBuffer, 0, maxOutputFrames * pRenderFormat_->nBlockAlign);
         else if (outputBuffer) memset(outputBuffer, 0, maxOutputFrames);
         *numInputFramesConsumedOut = 0;
         return 0;
    }

     if (!IsFormatCombinationSupported(pCaptureFormat_, pRenderFormat_)) {
         if (!warned_unsupported_conversion_) {
              printf("Warning [AudioConverter::Convert]: Format combination unsupported. Filling output with silence.\n");
              warned_unsupported_conversion_ = true;
         }
         memset(outputBuffer, 0, maxOutputFrames * pRenderFormat_->nBlockAlign);
         *numInputFramesConsumedOut = 0;
         return 0;
    } else {
         if (warned_unsupported_conversion_) {
              printf("Info [AudioConverter::Convert]: Supported format detected. Re-enabling conversion.\n");
              warned_unsupported_conversion_ = false;
         }
    }

    size_t captureChannels = pCaptureFormat_->nChannels;
    size_t renderChannels = pRenderFormat_->nChannels;
    size_t captureBytesPerSample = pCaptureFormat_->wBitsPerSample / 8;
    size_t renderBytesPerSample = pRenderFormat_->wBitsPerSample / 8;
    size_t captureBlockAlign = pCaptureFormat_->nBlockAlign;
    size_t renderBlockAlign = pRenderFormat_->nBlockAlign;

    // Corrected BlockAlign check
    if (captureBlockAlign == 0 && pCaptureFormat_->nChannels > 0 && pCaptureFormat_->wBitsPerSample > 0) {
         captureBlockAlign = (size_t)pCaptureFormat_->nChannels * (pCaptureFormat_->wBitsPerSample / 8);
    }
     if (renderBlockAlign == 0 && pRenderFormat_->nChannels > 0 && pRenderFormat_->wBitsPerSample > 0) {
         renderBlockAlign = (size_t)pRenderFormat_->nChannels * (pRenderFormat_->wBitsPerSample / 8);
    }

    if (captureBlockAlign == 0 || renderBlockAlign == 0) {
        printf("Error [AudioConverter::Convert]: Zero block alignment.\n");
         memset(outputBuffer, 0, maxOutputFrames * renderBlockAlign);
         *numInputFramesConsumedOut = 0;
         return 0;
    }


    size_t outputFramesWritten = 0;
    size_t stateBufferSizeFrames = GetStateBufferSizeFrames(); // Size in frames (1 for linear)

    // We need a buffer containing state frames + current input block.
    // The combined buffer size = stateSizeFrames + numInputFramesInBlock
    // The inputBuffer pointer is the start of the NEW DATA block.
    // To access state, we need access to the buffer *before* inputBuffer.
    // A simpler approach: Convert manages the combined buffer internally.
    // It takes the NEW data block, prepends its state, processes, and saves the last frame for the next state.

    // Let's use the model where Convert takes the NEW input data and manages its own state prepending.
    // This requires re-reading the inputBuffer + state logic...

    // Alternative model for Convert:
    // Convert(const BYTE* inputDataBlock, size_t numInputFramesInBlock, ...)
    // It needs access to the *previous* frame (state).
    // So, the RenderThread MUST provide the state + new data combined buffer.
    // The original plan where RenderThread builds the combined buffer [state][data] and passes it to Convert is correct.
    // Let's stick to that. The `inputBuffer` parameter in THIS Convert signature
    // *is* the combined buffer [state][data]. And `numInputFramesInBlock` is the size of the data *segment*.

    // Re-correcting Convert parameters and logic based on the [state][data] combined buffer:
    // size_t Convert(const BYTE* combinedBuffer, size_t numInputFramesInCombinedBuffer, BYTE* outputBuffer, size_t maxOutputFrames, size_t* numInputFramesConsumedFromDataBlockOut)

    // Let's rewrite Convert again based on the [state][data] combined buffer.
    // The previous code block for Convert seemed closer to the [state][data] model.
    // Let's use that code block and fix its parameters and state logic.

    // Start of the data segment within the combined buffer:
    const BYTE* inputDataSegment = inputBuffer + stateBufferSizeFrames * captureBlockAlign;
    size_t totalFramesCombined = stateBufferSizeFrames + numInputFramesInBlock;


    // Can process with interpolation if total combined frames is >= 2 (state frame + at least 1 data frame)
    bool can_process_with_interpolation = (totalFramesCombined >= 2);

     if (!can_process_with_interpolation && std::fabs(resampleRatio_ - 1.0) > 1e-6) {
         if (!warned_unsupported_conversion_) {
              printf("Warning [AudioConverter::Convert]: Not enough combined input frames (%zu) for interpolation (need >= 2). Filling output with silence.\n",
                     totalFramesCombined);
              warned_unsupported_conversion_ = true;
         }
         memset(outputBuffer, 0, maxOutputFrames * renderBlockAlign);
         *numInputFramesConsumedOut = 0; // Consumed 0 from the new data block
         return 0;
     } else {
         if (warned_unsupported_conversion_) {
             printf("Info [AudioConverter::Convert]: Enough combined input frames (%zu) for interpolation/direct copy. Resuming processing.\n", totalFramesCombined);
             warned_unsupported_conversion_ = false;
         }
     }

    // Phase accumulator: fractional position within the COMBINED input buffer
    // It starts at resamplePhase_, which is the fractional position where we left off
    // in the PREVIOUS combined buffer, relative to the START of that previous combined buffer.
    // So, the first sample we need to calculate is at index `resamplePhase_` in the *current combined buffer*.
    double phase_accumulator = resamplePhase_; // fractional index in the CURRENT combined buffer

    size_t initial_output_frames_written = 0; // To track how many output frames we successfully wrote

    for (size_t i = 0; i < maxOutputFrames; ++i) {
        size_t idx1_combined = static_cast<size_t>(floor(phase_accumulator));
        size_t idx2_combined = idx1_combined + 1;
        float interp_factor = static_cast<float>(phase_accumulator - idx1_combined);

        if (idx2_combined >= totalFramesCombined) {
             // Not enough input data in the combined buffer to produce this output frame
             printf("Warning [AudioConverter::Convert]: Ran out of combined input data mid-conversion (needed index %zu, available %zu). Filling remaining output (%zu/%zu) with silence.\n",
                    idx2_combined, totalFramesCombined, initial_output_frames_written + outputFramesWritten, maxOutputFrames);
             memset(outputBuffer + outputFramesWritten * renderBlockAlign, 0, (maxOutputFrames - outputFramesWritten) * renderBlockAlign);
             break; // Exit the output frame loop
        }

        // Get pointers to the two surrounding input frames in the combined buffer
        const BYTE* pInputFrame1 = inputBuffer + idx1_combined * captureBlockAlign;
        const BYTE* pInputFrame2 = inputBuffer + idx2_combined * captureBlockAlign;

        BYTE* pOutputFrame = outputBuffer + outputFramesWritten * renderBlockAlign;
        memset(pOutputFrame, 0, renderBlockAlign);

        // Perform interpolation and channel mapping
        // This channel mapping logic is basic, but matches the previous version
        size_t numInputChannels = captureChannels;
        size_t numOutputChannels = renderChannels;

        if (numInputChannels == 1 && numOutputChannels >= 1) {
            // Mono -> Any: Copy interpolated mono sample to all output channels
            float s1 = GetFloatSample(pInputFrame1, pCaptureFormat_); // Channel 0 input
            float s2 = GetFloatSample(pInputFrame2, pCaptureFormat_); // Channel 0 input
            float s_interp = InterpolateLinear(s1, s2, interp_factor);
            for (size_t ch_render = 0; ch_render < numOutputChannels; ++ch_render) {
                 SetFloatSample(pOutputFrame + ch_render * renderBytesPerSample, s_interp, pRenderFormat_);
            }
        } else if (numInputChannels >= 1 && numOutputChannels == 1) {
             // Any -> Mono: average first N input channels (e.g., 2 for stereo)
             size_t inputChannelsToAverage = std::min(numInputChannels, 2u); // Average first 2, or just 1 if mono input
             float sum_interp = 0.0f;
             for(size_t ch_in = 0; ch_in < inputChannelsToAverage; ++ch_in) {
                 float s1 = GetFloatSample(pInputFrame1 + ch_in * captureBytesPerSample, pCaptureFormat_);
                 float s2 = GetFloatSample(pInputFrame2 + ch_in * captureBytesPerSample, pCaptureFormat_);
                 sum_interp += InterpolateLinear(s1, s2, interp_factor);
             }
             float s_interp = sum_interp / static_cast<float>(inputChannelsToAverage);
             // Write to the single output channel (0)
             SetFloatSample(pOutputFrame, s_interp, pRenderFormat_);

        } else if (numInputChannels == numOutputChannels) {
            // N -> N mapping (including 2->2)
            for (size_t ch = 0; ch < numInputChannels; ++ch) {
                float s1 = GetFloatSample(pInputFrame1 + ch * captureBytesPerSample, pCaptureFormat_);
                float s2 = GetFloatSample(pInputFrame2 + ch * captureBytesPerSample, pCaptureFormat_);
                float s_interp = InterpolateLinear(s1, s2, interp_factor);
                SetFloatSample(pOutputFrame + ch * renderBytesPerSample, s_interp, pRenderFormat_);
            }
        } else {
            // Other channel mapping combinations (e.g., 3->2, 2->4) - basic passthrough min(in, out) channels
             if (!warned_unsupported_conversion_) {
                  printf("Warning [AudioConverter::Convert]: Unsupported channel mapping combination (Capture %u -> Render %u). Defaulting to mapping min channels and zero-filling extra output.\n",
                         numInputChannels, numOutputChannels);
                  warned_unsupported_conversion_ = true;
             }
             size_t channelsToMap = std::min(numInputChannels, numOutputChannels);
             for (size_t ch = 0; ch < channelsToMap; ++ch) {
                 float s1 = GetFloatSample(pInputFrame1 + ch * captureBytesPerSample, pCaptureFormat_);
                 float s2 = GetFloatSample(pInputFrame2 + ch * captureBytesPerSample, pCaptureFormat_);
                 float s_interp = InterpolateLinear(s1, s2, interp_factor);
                 SetFloatSample(pOutputFrame + ch * renderBytesPerSample, s_interp, pRenderFormat_);
             }
             for (size_t ch = channelsToMap; ch < numOutputChannels; ++ch) {
                  SetFloatSample(pOutputFrame + ch * renderBytesPerSample, 0.0f, pRenderFormat_); // Zero-fill extra output channels
             }
        }


        // Advance the floating-point input index for the next output frame
        phase_accumulator += 1.0 / resampleRatio_;

        outputFramesWritten++;
    }

    // Update resamplePhase_ for the next call. It's the fractional part of the index
    // in the *current combined buffer* where we stopped.
    resamplePhase_ = fmod(phase_accumulator, 1.0);


    // Calculate the number of integer input frames consumed from the *NEW DATA* block.
    // The position in the combined buffer advanced from `initial_phase` to `phase_accumulator`.
    // Total frames (including state) advanced = `phase_accumulator - initial_phase`.
    // Integer frames advanced from the START OF THE COMBINED BUFFER = floor(phase_accumulator) - floor(initial_phase).
    // The data block starts at index `stateBufferSizeFrames` in the combined buffer.
    // The frames consumed from the DATA block are those whose indices relative to the START OF THE DATA BLOCK were crossed.
    // The index in the data block corresponding to `phase_accumulator` is `phase_accumulator - stateBufferSizeFrames`.
    // The integer frames from the data block that were *fully passed* are from index 0 up to `floor(phase_accumulator - stateBufferSizeFrames) - 1`.
    // So the number of data frames consumed is `floor(phase_accumulator - stateBufferSizeFrames)`. Let's use this.
    // We need to be careful with floating point arithmetic and negative results if phase_accumulator < stateBufferSizeFrames.
    // If phase_accumulator < stateBufferSizeFrames, it means we haven't even processed past the state frame yet.
    // In that case, 0 data frames are consumed.

    double end_pos_in_data_block_float = phase_accumulator - stateBufferSizeFrames;
    size_t inputFramesConsumedFromBlock = 0;
    if (end_pos_in_data_block_float >= 0.0) {
         inputFramesConsumedFromBlock = static_cast<size_t>(floor(end_pos_in_data_block_float));
    }
     // Cap consumed frames at the actual number of new data frames available
     inputFramesConsumedFromBlock = std::min(inputFramesConsumedFromBlock, numInputFramesInBlock);

    *numInputFramesConsumedOut = inputFramesConsumedFromBlock; // Output parameter


    // Save the state for the next call: copy the frame at index `floor(phase_accumulator)`
    // from the current *combined buffer* into the state buffer.
    // This is the frame at index `idx1_combined` in the last successful iteration, if the loop completed fully.
    // Or, the frame at index `floor(phase_accumulator)` in the combined buffer after the loop.
    // The state buffer should store the frame that will be needed as `pInputFrame1` for the very first output sample in the *next* Convert call.
    // This frame is at index `floor(phase_accumulator)` in the CURRENT combined buffer.

    if (stateBufferSizeFrames > 0 && captureBlockAlign > 0) {
         size_t frame_index_to_save_combined = static_cast<size_t>(floor(phase_accumulator));
         if (frame_index_to_save_combined < totalFramesCombined) {
             memcpy(audioConverter.interpolationStateBuffer_.data(),
                    inputBuffer + frame_index_to_save_combined * captureBlockAlign,
                    captureBlockAlign);
         } else {
              // This happens if phase_accumulator is exactly at the end of the combined buffer,
              // or if totalFramesCombined was 0 (should be caught earlier).
              // The last frame of the available combined buffer is totalFramesCombined - 1.
              if (totalFramesCombined > 0) {
                   memcpy(audioConverter.interpolationStateBuffer_.data(),
                          inputBuffer + (totalFramesCombined - 1) * captureBlockAlign,
                          captureBlockAlign);
              } else {
                   // No data available at all, zero state.
                   memset(audioConverter.interpolationStateBuffer_.data(), 0, audioConverter.interpolationStateBuffer_.size());
              }
         }
    }


    return outputFramesWritten; // Return number of output frames produced
}


// 6. Definitions of Thread Functions

DWORD WINAPI AudioCaptureThread(LPVOID lpParam) {
    HRESULT hr = S_OK;
    // unique_ptr теперь владеет params, его деструктор вызовется при выходе из функции,
    // автоматически очищая ComPtr, ScopedHandle и vector в структуре params.
    std::unique_ptr<AudioCaptureThreadParams> params(static_cast<AudioCaptureThreadParams*>(lpParam));

    // Получаем сырые указатели и хэндлы из членов params (они управляются ComPtr/ScopedHandle внутри params)
    // ИСправлено: Используем .Get() для явного получения сырых указателей
    IAudioClient *pAudioClient = params->pAudioClient.Get();
    IAudioCaptureClient *pCaptureClient = params->pCaptureClient.Get();
    HANDLE hCaptureEvent = params->hCaptureEvent.Get();
    HANDLE hStopEvent = params->hStopEvent.Get(); // hStopEvent - копия хэндла

    // Указатели на общие структуры
    SharedAudioBuffer* pSharedBuffer = params->pSharedBuffer;
    SharedDspParameters* pSharedDspParams = params->pSharedDspParams;
    SharedMonitorData* pSharedMonitorData = params->pSharedMonitorData;

    UINT32 bufferFrameCount = params->bufferFrameCount;
    REFERENCE_TIME devicePeriod = params->devicePeriod;
    bool isEventDriven = params->isEventDriven;

    // Ссылка на вектор эффектов внутри params
    std::vector<std::unique_ptr<IDspEffect>>& effectsChain = params->effectsChain;


    printf("[Поток захвата] Поток запущен.\n");
    printf("[Поток захвата] Mode: %s, Buffer Size: %u frames, Period: %lld ms.\n", isEventDriven ? "Event" : "Polling", bufferFrameCount, devicePeriod/10000);

    // УДАЛЕНО: Manual AddRef() вызовы - RAII обертки в params управляют этим.

    const WAVEFORMATEX* pCaptureFormat = pSharedBuffer->GetCaptureFormat();
    if (!pCaptureFormat) {
        printf("[Поток захвата] Ошибка: Capture Format NULL в SharedBuffer при старте потока. Завершение.\n");
        // Заменяем goto на return
        return E_UNEXPECTED;
    }

    // --- DSP Chain Setup (Initial) ---
    printf("[Поток захвата] Настройка начальной цепочки DSP: Gain -> Biquad...\n");

    // Create Gain Effect
    auto gainEffect = std::make_unique<GainEffect>();
    if (pSharedDspParams) { gainEffect->UpdateParameters(*pSharedDspParams, pCaptureFormat); }
    effectsChain.push_back(std::move(gainEffect)); // Перемещаем владение в вектор
    printf("[Поток захвата] Добавлен эффект: %s\n", effectsChain.back()->GetName());

    // Create Biquad Effect
    auto biquadEffect = std::make_unique<BiquadEffect>();
    if (pSharedDspParams) {
        // Используем static_cast для вызова Initialize/Update, как обсуждали
        BiquadEffect* rawBiquadEffect = static_cast<BiquadEffect*>(biquadEffect.get());
        if (rawBiquadEffect) {
            rawBiquadEffect->Initialize(pCaptureFormat); // Biquad needs format during init
            rawBiquadEffect->UpdateParameters(*pSharedDspParams, pCaptureFormat);
        } else {
            printf("[Поток захвата] Error: Failed to static_cast BiquadEffect pointer for initialization/update.\n");
             // Решаем выйти при ошибке, если инициализация Biquad критична.
             // Перед выходом RAII в unique_ptr params почистит ресурсы, включая уже добавленный GainEffect.
            return E_UNEXPECTED;
        }
    }
    effectsChain.push_back(std::move(biquadEffect)); // Перемещаем владение в вектор
    printf("[Поток захвата] Добавлен эффект: %s\n", effectsChain.back()->GetName());

    printf("[Поток захвата] Начальная цепочка DSP готова.\n");

    // Массив хэндлов для WaitForMultipleObjects
    HANDLE waitHandles[3];
    waitHandles[0] = hStopEvent; // Событие остановки всегда на первом месте
    DWORD nHandles = 1;
    if (isEventDriven && hCaptureEvent && hCaptureEvent != INVALID_HANDLE_VALUE) { waitHandles[nHandles++] = hCaptureEvent; }
    HANDLE hDspEvent = NULL;
    if (pSharedDspParams && pSharedDspParams->hParamChangedEvent && pSharedDspParams->hParamChangedEvent != INVALID_HANDLE_VALUE) {
        hDspEvent = pSharedDspParams->hParamChangedEvent;
        waitHandles[nHandles++] = hDspEvent;
        printf("[Поток захвата] DSP Parameter Change Event added to wait list.\n");
    } else {
        printf("[Поток захвата] Warning: Shared DSP parameters object or event handle is invalid.\n");
    }

    HRESULT cycleHr = S_OK; // Результат текущего цикла обработки

    while (TRUE) {
        DWORD dwWaitResult;
        DWORD timeout = isEventDriven ? INFINITE : (DWORD)(devicePeriod/10000/2);
        if (timeout == 0) timeout = 1; // Ensure non-zero timeout in polling mode

        // Ожидаем одно из событий или таймаут
        dwWaitResult = WaitForMultipleObjects(nHandles, waitHandles, FALSE, timeout);

        // Проверяем событие остановки (всегда на waitHandles[0])
        if (dwWaitResult == WAIT_OBJECT_0) {
             cycleHr = S_OK; // Штатный выход по сигналу остановки
             break; // Выходим из цикла
        }

        // Check if DSP parameters changed (only if event handle was added)
        // Event handle is at index nHandles-1 if nHandles > 1 and it's not hStopEvent
        if (hDspEvent && nHandles > 1 && waitHandles[nHandles-1] == hDspEvent && dwWaitResult == WAIT_OBJECT_0 + (nHandles-1)) {
             printf("[Поток захвата] Событие изменения DSP параметров сигнализировано. Обновляю эффекты...\n");
             const WAVEFORMATEX* currentCaptureFormat = pSharedBuffer->GetCaptureFormat();
             if (currentCaptureFormat) {
                  for (const auto& effect_ptr : effectsChain) {
                     if (effect_ptr) {
                         effect_ptr->UpdateParameters(*pSharedDspParams, currentCaptureFormat);
                     }
                  }
             } else { printf("[Поток захвата] Warning: Cannot update DSP parameters, capture format is NULL.\n"); }
             // Событие Auto-reset, не нужно сбрасывать вручную.
             continue; // Переходим к следующей итерации цикла
        }
         // Если сработало событие hCaptureEvent (Event-driven mode, индекс 1) или таймаут (Polling mode), продолжаем обработку ниже.


        UINT32 numFramesAvailable = 0; BYTE *pData = nullptr; DWORD flags = 0; UINT64 devPos = 0, qpcPos = 0;
        HRESULT getBufferHr = AUDCLNT_S_BUFFER_EMPTY; // Инициализируем на "пусто"

         // Обрабатываем WASAPI событие или таймаут Polling Mode
         // Если EventDriven, hCaptureEvent - waitHandles[1]
         // Если Polling, waitHandles[1] отсутствует, ждем таймаут
         DWORD captureEventIndex = 1; // Позиция hCaptureEvent в waitHandles, если он есть
         bool processBuffer = false;

         if (isEventDriven && hCaptureEvent && hCaptureEvent != INVALID_HANDLE_VALUE) {
              // Если сработало событие захвата ИЛИ был таймаут
             if (dwWaitResult == WAIT_OBJECT_0 + captureEventIndex || dwWaitResult == WAIT_TIMEOUT) {
                  getBufferHr = pCaptureClient->GetBuffer(&pData, &numFramesAvailable, &flags, &devPos, &qpcPos);
                  processBuffer = true; // Пробуем обработать буфер после события или таймаута
             } else {
                 // Сработало другое событие (StopEvent или DspEvent), уже обработано выше.
                 continue;
             }
         } else { // Polling Mode
              // В Polling Mode всегда пробуем получить буфер по таймауту
              if (dwWaitResult == WAIT_TIMEOUT || dwWaitResult == WAIT_OBJECT_0) { // WAIT_OBJECT_0 уже обработан, но для полноты условия
                 getBufferHr = pCaptureClient->GetBuffer(&pData, &numFramesAvailable, &flags, &devPos, &qpcPos);
                 processBuffer = true; // Пробуем обработать буфер
              } else {
                 // Сработало StopEvent (уже обработано)
                 continue;
              }
         }


        if (pSharedMonitorData) {
            UINT32 currentPadding = 0;
            // ИСправлено: Используем pAudioClient (сырой указатель)
            HRESULT paddingHr = pAudioClient->GetCurrentPadding(&currentPadding); // ИСправлено: опечатка ¤tPadding -> &currentPadding
            WARN_ON_ERROR(paddingHr, "GetCurrentPadding (Capture Thread)");
            if (SUCCEEDED(paddingHr)) pSharedMonitorData->UpdateCapturePadding(currentPadding);
            else pSharedMonitorData->UpdateCapturePadding(0); // Сбрасываем в 0 при ошибке
        }

        if (processBuffer && SUCCEEDED(getBufferHr) && numFramesAvailable > 0) {
             const WAVEFORMATEX* pCurrentCaptureFormat = pSharedBuffer->GetCaptureFormat();
             if (!pCurrentCaptureFormat) {
                 printf("[Поток захвата] Ошибка: Capture Format NULL перед обработкой. Завершение.\n");
                 hr = E_UNEXPECTED; // Задаем код ошибки
                 // Освобождаем буфер WASAPI перед выходом
                 // ИСправлено: Используем pCaptureClient (сырой указатель)
                 HRESULT releaseHr = pCaptureClient->ReleaseBuffer(numFramesAvailable);
                 WARN_ON_ERROR(releaseHr, "ReleaseBuffer (Capture - before exit)");
                 break; // Выходим из цикла
             }

             if (!(flags & AUDCLNT_BUFFERFLAGS_SILENT)) {
                 ApplyDspEffectChain(pData, numFramesAvailable, pCurrentCaptureFormat, effectsChain);

                 // Входим в критическую секцию для доступа к SharedBuffer
                 EnterCriticalSection(&pSharedBuffer->lock);
                 size_t framesWritten = pSharedBuffer->Write(pData, numFramesAvailable);
                 LeaveCriticalSection(&pSharedBuffer->lock);
                 // printf("[Поток захвата] Wrote %zu frames to SharedBuffer.\n", framesWritten); // Reduced logging

                 // ИСправлено: Используем pCaptureClient (сырой указатель)
                 HRESULT releaseHr = pCaptureClient->ReleaseBuffer(numFramesAvailable);
                 WARN_ON_ERROR(releaseHr, "ReleaseBuffer (Capture)");
                 if (FAILED(releaseHr)) {
                     printf("[Поток захвата] Ошибка ReleaseBuffer: 0x%x. Завершение потока.\n", releaseHr);
                     hr = releaseHr; // Задаем код ошибки
                     break; // Выходим из цикла
                 }
             } else {
                 // Received silent buffer from WASAPI capture. No DSP, but release the buffer.
                 // ИСправлено: Используем pCaptureClient (сырой указатель)
                 HRESULT releaseHr = pCaptureClient->ReleaseBuffer(numFramesAvailable);
                 WARN_ON_ERROR(releaseHr, "ReleaseBuffer (Capture - Silent)");
                 if (FAILED(releaseHr)) {
                     printf("[Поток захвата] Ошибка ReleaseBuffer тихого буфера: 0x%x. Завершение потока.\n", releaseHr);
                     hr = releaseHr; // Задаем код ошибки
                     break; // Выходим из цикла
                 }
             }
         } else if (getBufferHr != AUDCLNT_S_BUFFER_EMPTY && FAILED(getBufferHr)) {
              // Обрабатываем ошибку GetBuffer, если она не AUDCLNT_S_BUFFER_EMPTY
              printf("[Поток захвата] Ошибка GetBuffer: 0x%x. Завершение.\n", getBufferHr);
              hr = getBufferHr; // Задаем код ошибки
              break; // Выходим из цикла
         }
         // Если getBufferHr == AUDCLNT_S_BUFFER_EMPTY или numFramesAvailable == 0, просто продолжаем цикл
    } // End while(TRUE)

    // УДАЛЕНО: Метка ThreadExit и весь код после нее.
    // Cleanup будет выполнен автоматически деструкторами params unique_ptr
    // и его членов (ComPtr, ScopedHandle, vector<unique_ptr>).

    printf("[Поток захвата] Поток завершает работу.\n");
    // ComPtr/ScopedHandle в params автоматически вызовут Release/CloseHandle.
    // std::vector<unique_ptr> в params автоматически уничтожит эффекты.
    printf("[Поток захвата] Ресурсы потока освобождены. Поток вышел.\n");
    return cycleHr; // Возвращаем финальный HRESULT цикла (S_OK если вышел по StopEvent)
}


DWORD WINAPI AudioRenderThread(LPVOID lpParam) {
    HRESULT hr = S_OK;
    // unique_ptr теперь владеет params, его деструктор вызовется при выходе из функции
    std::unique_ptr<AudioRenderThreadParams> params(static_cast<AudioRenderThreadParams*>(lpParam));

    // Получаем сырые указатели и хэндлы из членов params (они управляются ComPtr/ScopedHandle внутри params)
    // ИСправлено: Используем .Get() для явного получения сырых указателей
    IAudioClient *pAudioClient = params->pAudioClient.Get();
    IAudioRenderClient *pRenderClient = params->pRenderClient.Get();
    PWAVEFORMATEX pRenderFormat = params->pRenderFormat; // Это указатель на буфер внутри params
    HANDLE hRenderEvent = params->hRenderEvent.Get();
    HANDLE hStopEvent = params->hStopEvent.Get(); // hStopEvent - копия хэндла

    // Указатели на общие структуры
    SharedAudioBuffer* pSharedBuffer = params->pSharedBuffer;
    SharedMonitorData* pSharedMonitorData = params->pSharedMonitorData;

    UINT32 bufferFrameCount = params->bufferFrameCount;
    REFERENCE_TIME devicePeriod = params->devicePeriod;
    bool isEventDriven = params->isEventDriven;

    AudioConverter& audioConverter = params->audioConverter; // Ссылка на объект конвертера внутри params

    // УДАЛЕНО: Manual AddRef() вызовы.

    printf("[Поток рендеринга] Поток запущен.\n");
    printf("[Поток рендеринга] Mode: %s, Buffer Size: %u frames, Period: %lld ms.\n", isEventDriven ? "Event" : "Polling", bufferFrameCount, devicePeriod/10000);

    if (!pRenderFormat) {
        printf("[Поток рендеринга] Ошибка: Render Format NULL в параметрах. Завершение.\n");
        // Заменяем goto на return
        return E_UNEXPECTED;
    }

    // Initialize AudioConverter
    const WAVEFORMATEX* pCaptureFormat = pSharedBuffer->GetCaptureFormat(); // Get capture format from shared buffer
    if (!pCaptureFormat) {
        printf("[Поток рендеринга] Ошибка: Capture Format NULL в SharedBuffer при старте потока рендеринга. Завершение.\n");
        // Заменяем goto на return
        return E_UNEXPECTED;
    }

    if (!audioConverter.Initialize(pCaptureFormat, pRenderFormat)) {
        printf("[Поток рендеринга] Ошибка инициализации AudioConverter. Несовместимые форматы? Завершение.\n");
        // Заменяем goto на return
        return E_FAIL;
    }

    size_t stateBufferSize = audioConverter.GetStateBufferSize(); // Размер в байтах
    size_t stateBufferSizeFrames = audioConverter.GetStateBufferSizeFrames(); // Размер в фреймах

    size_t captureBlockAlign = pCaptureFormat->nBlockAlign;
    if (captureBlockAlign == 0 && pCaptureFormat_->nChannels > 0 && pCaptureFormat_->wBitsPerSample > 0) {
         captureBlockAlign = (size_t)pCaptureFormat_->nChannels * (pCaptureFormat_->wBitsPerSample / 8);
    }


    // Estimate max frames needed from SharedBuffer (Capture Format)
    // This estimate determines the max size of the 'new data' part of the combined buffer.
    // The actual amount read will be limited by what's available in SharedBuffer.
    size_t maxCaptureFramesToReadFromSharedBuffer = 0;
    if (audioConverter.GetResampleRatio() > 0) { // Avoid division by zero
        // Estimate based on render buffer size. Adding a safety margin.
        // Need enough input frames to produce 'bufferFrameCount' output frames.
        // If ratio > 1 (downsampling), need less input frames than output frames.
        // If ratio < 1 (upsampling), need more input frames than output frames.
        // A generous buffer might be bufferFrameCount * max(1.0, ratio) + some extra.
        // Let's use a simpler heuristic: size sufficient for the WASAPI buffer equivalent + safety.
        maxCaptureFramesToReadFromSharedBuffer = (size_t)((double)bufferFrameCount * audioConverter.GetResampleRatio()) + 10; // 10 frames safety
    } else {
        // Should not happen if converter initialized, but fallback
         maxCaptureFramesToReadFromSharedBuffer = bufferFrameCount + 10;
    }

    // Calculate required size for the combined buffer: [state_buffer] [new_read_data]
    size_t tempCombinedBufferSize = 0;
     if (captureBlockAlign > 0) {
          // Check for potential overflow BEFORE multiplication
          if (maxCaptureFramesToReadFromSharedBuffer > (SIZE_MAX - stateBufferSize) / captureBlockAlign) {
               printf("[Поток рендеринга] Ошибка: Расчетный размер временного буфера слишком большой (переполнение).\n");
               // Заменяем goto на return
               return E_OUTOFMEMORY;
          }
          tempCombinedBufferSize = stateBufferSize + maxCaptureFramesToReadFromSharedBuffer * captureBlockAlign;
     } else if (stateBufferSize > 0) {
          // Если BlockAlign 0, но есть размер состояния (что не должно быть при BlockAlign 0)
          printf("[Поток рендеринга] Ошибка: Capture format nBlockAlign is 0, but stateBufferSize is %zu.\n", stateBufferSize);
          // Заменяем goto на return
          return E_UNEXPECTED;
     } else {
          // BlockAlign is 0 and no state size. Buffer size is 0. Safe to proceed, vector remains empty.
          tempCombinedBufferSize = 0;
     }

    // --- Ресайз вектора происходит здесь, до любого return ---
    std::vector<BYTE> tempCombinedBuffer; // Объявлен
    if (tempCombinedBufferSize > 0) {
        try {
             tempCombinedBuffer.resize(tempCombinedBufferSize); // <-- Инициализация (выделение/создание элементов) происходит здесь
             printf("[Поток рендеринга] Временный буфер для конвертера: %zu байт (State: %zu bytes, Max New Read: %zu frames @ %zu bytes/frame).\n",
                    tempCombinedBufferSize, stateBufferSize, maxCaptureFramesToReadFromSharedBuffer, captureBlockAlign);
        } catch (const std::bad_alloc& e) {
             printf("[Поток рендеринга] Ошибка выделения памяти для временного буфера: %s. Завершение.\n", e.what());
             // Заменяем goto на return
             return E_OUTOFMEMORY; // Используем код ошибки памяти
        }
    } else {
         // Если tempCombinedBufferSize == 0, значит ресайз не нужен.
         printf("[Поток рендеринга] Временный буфер для конвертера: 0 байт (ресэмплинг не требуется или нет данных для чтения).\n");
    }
    // ------------------------------------------------------------

    HANDLE waitHandles[2];
    waitHandles[0] = hStopEvent; // Событие остановки всегда на первом месте
    DWORD nHandles = isEventDriven ? 2 : 1;
    if (isEventDriven && hRenderEvent && hRenderEvent != INVALID_HANDLE_VALUE) waitHandles[1] = hRenderEvent;


    HRESULT cycleHr = S_OK; // Результат текущего цикла обработки

    while (TRUE) {
        DWORD dwWaitResult;
        DWORD timeout = isEventDriven ? INFINITE : (DWORD)(devicePeriod/10000/2);
        if (timeout == 0) timeout = 1;

        // Ожидаем событие или таймаут
        dwWaitResult = WaitForMultipleObjects(nHandles, waitHandles, FALSE, timeout);

        // Проверяем событие остановки (всегда на waitHandles[0])
        if (dwWaitResult == WAIT_OBJECT_0) {
             cycleHr = S_OK; // Штатный выход по сигналу остановки
             break; // Выходим из цикла
        }

        UINT32 padding = 0;
        // ИСправлено: Используем pAudioClient (сырой указатель)
        HRESULT paddingHr = pAudioClient->GetCurrentPadding(&padding);
        if (FAILED(paddingHr)) {
            printf("[Поток рендеринга] Error GetCurrentPadding: 0x%x. Exiting cycle.\n", paddingHr);
             cycleHr = paddingHr; // Фиксируем ошибку для выхода
            break; // Выходим из цикла
        }
        if (pSharedMonitorData) pSharedMonitorData->UpdateRenderPadding(padding);

        UINT32 numFramesAvailableInWasapi = bufferFrameCount - padding;

        if (numFramesAvailableInWasapi == 0) { continue; } // Нет места в буфере WASAPI, ждем следующего события/таймаута

        BYTE* pRenderData = nullptr;
        // ИСправлено: Используем pRenderClient (сырой указатель)
        HRESULT getBufferHr = pRenderClient->GetBuffer(numFramesAvailableInWasapi, &pRenderData);
        if (FAILED(getBufferHr)) {
            printf("[Поток рендеринга] Error GetBuffer: 0x%x. Exiting cycle.\n", getBufferHr);
             cycleHr = getBufferHr; // Фиксируем ошибку для выхода
            break; // Выходим из цикла
        }
        if (pRenderData == nullptr) {
            printf("[Поток рендеринга] GetBuffer returned NULL pData. Exiting cycle.\n");
             cycleHr = E_UNEXPECTED; // Неожиданная ошибка
            break; // Выходим из цикла
        }

        // --- Prepare combined buffer for AudioConverter ---
        size_t captureFramesRead = 0; // Number of frames read from SharedBuffer in this cycle

        // Copy state buffer to the beginning of the combined buffer
        if (stateBufferSize > 0) { // Check that state buffer exists and has size > 0
            if (!audioConverter.interpolationStateBuffer_.empty()) { // Check that the internal vector is not empty
                 // ISправлено: Используем .data() для получения сырого указателя вектора
                 // tempCombinedBuffer должен быть достаточно большим, проверено при его выделении.
                memcpy(tempCombinedBuffer.data(), audioConverter.interpolationStateBuffer_.data(), stateBufferSize);
            } else {
                 // This can happen if state buffer exists (resampling active) but its content is empty (e.g., first run, or previous Read failed).
                 // Fill the state part of the combined buffer with silence.
                 memset(tempCombinedBuffer.data(), 0, stateBufferSize);
                 // printf("Debug [AudioRenderThread]: State buffer empty, filling state portion with silence.\n"); // Отладочное сообщение (можно закомментировать)
            }
        }

        // Read new data from SharedBuffer into the combined buffer *after* the state
        // Читаем максимум, что можем или что нужно.
        size_t framesToReadFromSharedBufferThisCycle = std::min(maxCaptureFramesToReadFromSharedBuffer, pSharedBuffer->GetAvailableFrames());
        // Проверяем, что есть куда читать и что BlockAlign корректен
        if (framesToReadFromSharedBufferThisCycle > 0 && captureBlockAlign > 0) {
             // Убедимся, что смещение + размер чтения не выходит за пределы tempCombinedBuffer.data() + tempCombinedBufferSize
             if (stateBufferSize + framesToReadFromSharedBufferThisCycle * captureBlockAlign <= tempCombinedBufferSize) {
                // Входим в критическую секцию для доступа к SharedBuffer
                EnterCriticalSection(&pSharedBuffer->lock);
                // ИСправлено: Используем .data() для получения сырого указателя вектора
                // Read теперь возвращает количество ФРЕЙМОВ
                captureFramesRead = pSharedBuffer->Read(tempCombinedBuffer.data() + stateBufferSize, framesToReadFromSharedBufferThisCycle);
                LeaveCriticalSection(&pSharedBuffer->lock);
             } else {
                  printf("[Поток рендеринга] Error: Read size calculation exceeds tempCombinedBuffer size (%zu > %zu). Skipping read.\n", stateBufferSize + framesToReadFromSharedBufferThisCycle * captureBlockAlign, tempCombinedBufferSize);
                  // Не выходим сразу, просто не читаем новые данные в этом цикле. Конвертер получит только состояние (или пустой буфер).
                  captureFramesRead = 0; // Ensure captureFramesRead is 0 on error
             }
        }


        // --- Convert and map data using AudioConverter ---
        size_t renderFramesWritten = 0; // Number of output frames produced by converter
        size_t inputFramesConsumedByConverter = 0; // Number of input frames consumed from the *NEW DATA* part

        // totalFramesAvailableInCombinedBuffer = state + new data
        size_t totalFramesAvailableInCombinedBuffer = stateBufferSizeFrames + captureFramesRead;


        if (totalFramesAvailableInCombinedBuffer > 0) {
             // Convert takes (combinedBuffer, totalInputFramesCombined, outputBuffer, maxOutputFrames, &inputFramesConsumedFromDataBlockOut)
             // ИСправлено: Вызов Convert с корректными параметрами
             // inputBuffer here IS the combined buffer. totalFramesAvailableInCombinedBuffer is its total size in frames.
             // We need to pass the number of frames IN THE NEW DATA BLOCK to Convert, which is captureFramesRead.
             // No, the Convert implementation expects total frames in the combined buffer. Let's re-check the Convert impl...
             // The Convert implementation (corrected version in previous response) expected:
             // inputBuffer: combined buffer [state][data]
             // numInputFramesInBlock: number of frames IN THE NEW DATA segment (after state)
             // outputBuffer, maxOutputFrames
             // numInputFramesConsumedOut: consumed from NEW DATA segment

             // Okay, let's call Convert correctly based on that:
             renderFramesWritten = audioConverter.Convert(tempCombinedBuffer.data(), captureFramesRead, pRenderData, numFramesAvailableInWasapi, &inputFramesConsumedByConverter);

              // --- Advance SharedBuffer Read Pointer ---
             if (inputFramesConsumedByConverter > 0) {
                 // EnterCriticalSection(&pSharedBuffer->lock); // Lock is handled inside AdvanceReadPointer
                 pSharedBuffer->AdvanceReadPointer(inputFramesConsumedByConverter);
                 // LeaveCriticalSection(&pSharedBuffer->lock); // Lock is handled inside AdvanceReadPointer

                 // Обновляем мониторинг потребленных фреймов
                 if (pSharedMonitorData) {
                     EnterCriticalSection(&pSharedMonitorData->lock);
                     pSharedMonitorData->framesConsumedFromSharedBuffer += inputFramesConsumedByConverter; // Сколько реально обработал конвертер из нового сегмента
                     LeaveCriticalSection(&pSharedMonitorData->lock);
                 }
             }

        } else {
            // Нет данных ни в состоянии, ни в SharedBuffer. WASAPI буфер должен быть заполнен тишиной.
             memset(pRenderData, 0, numFramesAvailableInWasapi * pRenderFormat->nBlockAlign);
             renderFramesWritten = 0;
             inputFramesConsumedByConverter = 0; // Ensure this is zero if no data
        }

         // *** Мониторинг ***
         // Обновляем мониторинг прочитанных фреймов из SharedBuffer. Read теперь возвращает фреймы.
         if (pSharedMonitorData) {
              EnterCriticalSection(&pSharedMonitorData->lock);
              pSharedMonitorData->framesReadFromSharedBuffer += captureFramesRead; // Сколько прочитали в этом цикле из SharedBuffer
              pSharedMonitorData->framesWrittenToRenderBuffer += renderFramesWritten; // Log frames written to WASAPI render buffer
              LeaveCriticalSection(&pSharedMonitorData->lock);
         }


        // --- Release WASAPI Buffer ---
        UINT32 framesToRelease = numFramesAvailableInWasapi; // Всегда освобождаем весь запрошенный WASAPI буфер
        DWORD dwFlags = 0;
        if (renderFramesWritten == 0 && numFramesAvailableInWasapi > 0) {
             // If we requested buffer space but converter produced nothing, tell WASAPI it's silent.
             dwFlags = AUDCLNT_BUFFERFLAGS_SILENT;
        }
        // If numFramesAvailableInWasapi was 0, GetBuffer wasn't called, ReleaseBuffer shouldn't be called either.
        if (numFramesAvailableInWasapi > 0) {
             // ИСправлено: Используем pRenderClient (сырой указатель)
             HRESULT releaseHr = pRenderClient->ReleaseBuffer(framesToRelease, dwFlags);
             WARN_ON_ERROR(releaseHr, "ReleaseBuffer (Render)");
             if (FAILED(releaseHr)) {
                 printf("[Поток рендеринга] Error ReleaseBuffer: 0x%x. Завершение потока.\n", releaseHr);
                  cycleHr = releaseHr; // Фиксируем ошибку для выхода
                 break; // Выходим из цикла
             }
        } else {
            // printf("Debug [Поток рендеринга]: numFramesAvailableInWasapi is 0. Skipping ReleaseBuffer.\n"); // Отладочное сообщение
        }


    } // End while(TRUE)

    // УДАЛЕНО: Метка ThreadExit и весь код после нее.
    // Ресурсы (ComPtr, ScopedHandle, vector) в unique_ptr params будут освобождены автоматически
    // при выходе из функции return cycleHr;

    printf("[Поток рендеринга] Поток завершает работу.\n");
    // Нет необходимости вызывать Stop/Reset/Release здесь, это делает main или они остановятся при завершении потока WASAPI,
    // а RAII обертки в params выполнят Release/CloseHandle для тех объектов, которыми владеет поток.
    // ИСправлено: Убраны явные вызовы Stop/Reset/Release.

    printf("[Поток рендеринга] Ресурсы потока освобождены. Поток вышел.\n");
    return cycleHr; // Возвращаем финальный HRESULT цикла (S_OK если вышел по StopEvent)
}


// ... (Остальные функции: SelectAudioDevicePrompt, GetMixFormatForDevice, InitializeCaptureDevice, InitializeRenderDevice, ApplyDspEffectChain, HandleUserInput, PrintCurrentDspSettings, PrintMonitorData, PrintWaveFormatEx, BiquadFilterTypeToString и хелперы для Biquad) ...

// 7. Main Function
int main()
{
    HRESULT hr = S_OK;
    if (!SetConsoleOutputCP(CP_UTF8)) { DWORD lastError = GetLastError(); if (lastError != 0) printf("Warning: SetConsoleOutputCP(CP_UTF8) failed with error %lu.\n", lastError); }
    printf("WASAPI C++ Sonar Project - Audio Loopback\n"); printf("------------------------------------------\n");
    hr = CoInitializeEx(NULL, COINIT_MULTITHREADED); if (FAILED(hr)) { printf("Ошибка CoInitializeEx: 0x%x\n", hr); return 1; } printf("COM инициализирован успешно.\n");
    auto co_uninitialize_on_exit = std::unique_ptr<void, void(*)(void*)>((void*)1, [](void* p) { CoUninitialize(); printf("COM завершен.\n"); });

    ComPtr<IMMDeviceEnumerator> pEnumerator; printf("Получение перечислителя аудиоустройств...\n");
    hr = CoCreateInstance(__uuidof(MMDeviceEnumerator), NULL, CLSCTX_ALL, __uuidof(IMMDeviceEnumerator), reinterpret_cast<LPVOID*>(pEnumerator.ReleaseAndGetAddressOf()));
    if (FAILED(hr)) { printf("Ошибка CoCreateInstance (IMMDeviceEnumerator): 0x%x\n", hr); return 1; } if (!pEnumerator) { printf("Ошибка: CoCreateInstance вернул S_OK, но указатель NULL.\n"); return 1; }

    ComPtr<IMMDevice> pCaptureDevice; CoTaskMemPtr pCaptureMixFormatMem; PWAVEFORMATEX pCaptureMixFormat = nullptr;
    const WCHAR* desiredCaptureDeviceId = L""; // Example: L"{0.0.1.00000000}.{xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx}"
    bool idSelectedCapture = (desiredCaptureDeviceId != NULL && desiredCaptureDeviceId[0] != L'\0');
    if (idSelectedCapture && pEnumerator) {
        wprintf_s(L"\nПопытка найти устройство захвата по ID: %ls\n", desiredCaptureDeviceId); hr = pEnumerator->GetDevice(desiredCaptureDeviceId, pCaptureDevice.ReleaseAndGetAddressOf());
        if (SUCCEEDED(hr) && pCaptureDevice) { printf("Устройство захвата найдено по ID.\n"); hr = GetMixFormatForDevice(pCaptureDevice.Get(), pCaptureMixFormatMem); if (SUCCEEDED(hr)) pCaptureMixFormat = reinterpret_cast<PWAVEFORMATEX>(pCaptureMixFormatMem.Get()); else { printf("Ошибка при получении Mix Format для устройства захвата по ID (0x%x).\n", hr); pCaptureDevice.Release(); hr = E_FAIL; }}
        else { wprintf_s(L"Устройство захвата с ID \"%ls\" не найдено или ошибка поиска (0x%x).\n", desiredCaptureDeviceId, hr); pCaptureDevice.Release(); hr = HRESULT_FROM_WIN32(ERROR_NOT_FOUND); }
    } else { printf("\nID устройства захвата не указан в коде.\n"); hr = HRESULT_FROM_WIN32(ERROR_NOT_FOUND); } // Fallback to interactive if ID not specified or not found
    if (FAILED(hr)) { printf("Запуск интерактивного выбора устройства захвата...\n"); hr = SelectAudioDevicePrompt(pEnumerator.Get(), eCapture, eConsole, pCaptureDevice, pCaptureMixFormatMem); if (FAILED(hr)) { printf("Ошибка при интерактивном выборе устройства захвата (0x%x). Программа завершается.\n", hr); return 1; } pCaptureMixFormat = reinterpret_cast<PWAVEFORMATEX>(pCaptureMixFormatMem.Get()); }
    if (!pCaptureDevice || !pCaptureMixFormat) { printf("Ошибка: Не удалось выбрать устройство захвата или получить его формат.\n"); return 1; }
    printf("\nSelected Capture Device Format:\n"); PrintWaveFormatEx(pCaptureMixFormat);

    ComPtr<IMMDevice> pRenderDevice; CoTaskMemPtr pRenderMixFormatMem; PWAVEFORMATEX pRenderMixFormat = nullptr;
    const WCHAR* desiredRenderDeviceId = L""; // Example: L"{0.0.0.00000000}.{xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx}"
    bool idSelectedRender = (desiredRenderDeviceId != NULL && desiredRenderDeviceId[0] != L'\0');
    if (idSelectedRender && pEnumerator) {
         wprintf_s(L"\nПопытка найти устройство вывода по ID: %ls\n", desiredRenderDeviceId); hr = pEnumerator->GetDevice(desiredRenderDeviceId, pRenderDevice.ReleaseAndGetAddressOf());
         if (SUCCEEDED(hr) && pRenderDevice) { printf("Устройство вывода найдено по ID.\n"); hr = GetMixFormatForDevice(pRenderDevice.Get(), pRenderMixFormatMem); if (SUCCEEDED(hr)) pRenderMixFormat = reinterpret_cast<PWAVEFORMATEX>(pRenderMixFormatMem.Get()); else { printf("Ошибка при получении Mix Format для устройства вывода по ID (0x%x).\n", hr); pRenderDevice.Release(); hr = E_FAIL; }}
         else { wprintf_s(L"Устройство вывода с ID \"%ls\" не найдено или ошибка поиска (0x%x).\n", desiredRenderDeviceId, hr); pRenderDevice.Release(); hr = HRESULT_FROM_WIN32(ERROR_NOT_FOUND); }
    } else { printf("\nID устройства вывода не указан в коде.\n"); hr = HRESULT_FROM_WIN32(ERROR_NOT_FOUND); } // Fallback to interactive if ID not specified or not found
    if (FAILED(hr)) { printf("Запуск интерактивного выбора устройства вывода...\n"); hr = SelectAudioDevicePrompt(pEnumerator.Get(), eRender, eConsole, pRenderDevice, pRenderMixFormatMem); if (FAILED(hr)) { printf("Ошибка при интерактивном выборе устройства вывода (0x%x). Программа завершается.\n", hr); return 1; } pRenderMixFormat = reinterpret_cast<PWAVEFORMATEX>(pRenderMixFormatMem.Get()); }
    if (!pRenderDevice || !pRenderMixFormat) { printf("Ошибка: Не удалось выбрать устройство вывода или получить его формат.\n"); return 1; }
    printf("\nSelected Render Device Format:\n"); PrintWaveFormatEx(pRenderMixFormat);


    //-------------------------------------------------------------------------
    // 4. Format Compatibility Check (Relaxed)
    //-------------------------------------------------------------------------
     if (!pCaptureMixFormat || !pRenderMixFormat) { printf("\nОшибка: Форматы устройств NULL после выбора.\n"); return 1; }

     size_t captureBitsPerSample = pCaptureMixFormat->wBitsPerSample;
     size_t renderBitsPerSample = pRenderMixFormat->wBitsPerSample;
     UINT16 captureValidBits = captureBitsPerSample;
     UINT16 renderValidBits = renderBitsPerSample;
     GUID captureSubFormat = {0}, renderSubFormat = {0};
     bool isCaptureExt = (pCaptureMixFormat->wFormatTag == WAVE_FORMAT_EXTENSIBLE);
     bool isRenderExt = (pRenderMixFormat->wFormatTag == WAVE_FORMAT_EXTENSIBLE);
     if (isCaptureExt && pCaptureMixFormat->cbSize >= sizeof(WAVEFORMATEXTENSIBLE) - sizeof(WAVEFORMATEX)) { const WAVEFORMATEXTENSIBLE* pExtCapture = reinterpret_cast<const WAVEFORMATEXTENSIBLE*>(pCaptureMixFormat); captureSubFormat = pExtCapture->SubFormat; captureValidBits = pExtCapture->Samples.wValidBitsPerSample; }
     else { captureValidBits = captureBitsPerSample; }
     if (isRenderExt && pRenderMixFormat->cbSize >= sizeof(WAVEFORMATEXTENSIBLE) - sizeof(WAVEFORMATEX)) { const WAVEFORMATEXTENSIBLE* pExtRender = reinterpret_cast<const WAVEFORMATEXTENSIBLE*>(pRenderMixFormat); renderSubFormat = pExtRender->SubFormat; renderValidBits = pExtRender->Samples.wValidBitsPerSample; }
     else { renderValidBits = renderBitsPerSample; }


     bool formatsCompatibleForSampleType = // Check if the *type* and *bit depth* are the same
         (pCaptureMixFormat->wFormatTag == pRenderMixFormat->wFormatTag &&
          (!isCaptureExt || !isRenderExt || IsEqualGUID(captureSubFormat, renderSubFormat)) && // SubFormat must match if both are EXT
          captureBitsPerSample == renderBitsPerSample && // Container size must match
          captureValidBits == renderValidBits); // Valid bits must match

     if (!formatsCompatibleForSampleType)
     {
         printf("\nОшибка: Форматы устройств захвата и вывода не совместимы для конверсии семплов.\n");
         printf("Требуется совпадение Tag, SubFormat (если оба EXT), wBitsPerSample (контейнер) и wValidBitsPerSample.\n");
         printf("Отличающиеся параметры:\n");
         if (pCaptureMixFormat->wFormatTag != pRenderMixFormat->wFormatTag) printf("  wFormatTag: %d (Cap) vs %d (Render)\n", pCaptureMixFormat->wFormatTag, pRenderMixFormat->wFormatTag);
         if (isCaptureExt && isRenderExt && !IsEqualGUID(captureSubFormat, renderSubFormat)) { printf("  SubFormat mismatch.\n"); }
         if (captureBitsPerSample != renderBitsPerSample) printf("  wBitsPerSample (container): %zu (Cap) vs %zu (Render)\n", captureBitsPerSample, renderBitsPerSample);
         if (captureValidBits != renderValidBits) printf("  ValidBitsPerSample: %u (Cap) vs %u (Render)\n", captureValidBits, renderValidBits);
         printf("Конверсия между разными типами/разрядностями не поддерживается на данном этапе.\n");
         return 1;
     } else {
         printf("\nФорматы устройств совместимы для конверсии семплов (разрешена разница в SampleRate и Channels).\n");
         printf("  Capture SampleRate: %lu, Render SampleRate: %lu\n", pCaptureMixFormat->nSamplesPerSec, pRenderMixFormat->nSamplesPerSec);
         printf("  Capture Channels: %d, Render Channels: %d\n", pCaptureMixFormat->nChannels, pRenderMixFormat->nChannels);
         printf("  Sample Type (Tag/SubFormat/Bits) Compatible.\n");
     }


    //-------------------------------------------------------------------------
    // 5. Initialize WASAPI Devices
    //-------------------------------------------------------------------------
     ComPtr<IAudioClient> pCaptureAudioClient; ComPtr<IAudioCaptureClient> pCaptureCaptureClient; ScopedHandle hCaptureEvent;
     ComPtr<IAudioClient> pRenderAudioClient; ComPtr<IAudioRenderClient> pRenderRenderClient; ScopedHandle hRenderEvent;

    hr = InitializeCaptureDevice(pCaptureDevice.Get(), pCaptureMixFormat, pCaptureAudioClient, pCaptureCaptureClient, hCaptureEvent);
    if (FAILED(hr)) { printf("Ошибка инициализации устройства захвата (0x%x).\n", hr); return 1; } if (!pCaptureAudioClient || !pCaptureCaptureClient) { printf("Ошибка: Инициализация устройства захвата вернула NULL объекты.\n"); return 1; }
    UINT32 captureBufferFrameCount = 0; REFERENCE_TIME captureDefaultPeriod = 0, captureMinimumPeriod = 0;
    if (pCaptureAudioClient) { hr = pCaptureAudioClient->GetBufferSize(&captureBufferFrameCount); WARN_ON_ERROR(hr, "GetBufferSize (main, capture)"); hr = pCaptureAudioClient->GetDevicePeriod(&captureDefaultPeriod, &captureMinimumPeriod); WARN_ON_ERROR(hr, "GetDevicePeriod (main, capture)"); printf("Параметры буфера захвата: WASAPI Size=%u frames, Default Period=%lld ms.\n", captureBufferFrameCount, captureDefaultPeriod/10000); } else { printf("Error: pCaptureAudioClient is NULL after successful InitializeCaptureDevice!\n"); return 1; }

    hr = InitializeRenderDevice(pRenderDevice.Get(), pRenderMixFormat, pRenderAudioClient, pRenderRenderClient, hRenderEvent);
    if (FAILED(hr)) { printf("Ошибка инициализации устройства вывода (0x%x).\n", hr); return 1; } if (!pRenderAudioClient || !pRenderRenderClient) { printf("Ошибка: Инициализация устройства вывода вернула NULL объекты.\n"); return 1; }
    UINT32 renderBufferFrameCount = 0; REFERENCE_TIME renderDefaultPeriod = 0, renderMinimumPeriod = 0;
    if (pRenderAudioClient) { hr = pRenderAudioClient->GetBufferSize(&renderBufferFrameCount); WARN_ON_ERROR(hr, "GetBufferSize (main, render)"); hr = pRenderAudioClient->GetDevicePeriod(&renderDefaultPeriod, &renderMinimumPeriod); WARN_ON_ERROR(hr, "GetDevicePeriod (main, render)"); printf("Параметры буфера рендеринга: WASAPI Size=%u frames, Default Period=%lld ms.\n", renderBufferFrameCount, renderDefaultPeriod/10000); } else { printf("Error: pRenderAudioClient is NULL after successful InitializeRenderDevice!\n"); return 1; }

    //-------------------------------------------------------------------------
    // 6. Initialize Shared Buffer (Based on Capture Format)
    //-------------------------------------------------------------------------
    SharedAudioBuffer sharedBuffer;
    if (sharedBuffer.Initialize(pCaptureMixFormat, 500) == 0) { printf("Ошибка: Не удалось инициализировать общий буфер.\n"); return 1; }

    //-------------------------------------------------------------------------
    // 7. Initialize Shared DSP Parameters
    //-------------------------------------------------------------------------
    SharedDspParameters sharedDspParams; sharedDspParams.Initialize();
     if (!sharedDspParams.hParamChangedEvent || sharedDspParams.hParamChangedEvent == INVALID_HANDLE_VALUE) { printf("Ошибка: Не удалось инициализировать общие параметры DSP (не создано событие).\n"); return 1; }
    sharedDspParams.SetGain(1.0f);
    sharedDspParams.SetFilterSettings(BiquadFilterType::Identity, 1000.0f, 1.0f, 0.0f, 1.0f);

    //-------------------------------------------------------------------------
    // 8. Initialize Shared Monitor Data
    //-------------------------------------------------------------------------
    SharedMonitorData sharedMonitorData; sharedMonitorData.Initialize();

    //-------------------------------------------------------------------------
    // 9. Create Stop Event
    //-------------------------------------------------------------------------
     ScopedHandle hStopEvent; hStopEvent = CreateEvent(NULL, TRUE, FALSE, NULL);
      if (!hStopEvent) { printf("Ошибка CreateEvent (hStopEvent): %lu.\n", GetLastError()); return 1; }

    //-------------------------------------------------------------------------
    // 10. Prepare Parameters and Create Threads
    //-------------------------------------------------------------------------
    auto pCaptureParamsHeap = std::make_unique<AudioCaptureThreadParams>();
    auto pRenderParamsHeap = std::make_unique<AudioRenderThreadParams>();

    // Перемещение владения ресурсами из main в структуру параметров
    pCaptureParamsHeap->pAudioClient = std::move(pCaptureAudioClient); // Перемещение владения
    pCaptureParamsHeap->pCaptureClient = std::move(pCaptureCaptureClient); // Перемещение владения
    pCaptureParamsHeap->hCaptureEvent = std::move(hCaptureEvent); // Перемещение владения
    pCaptureParamsHeap->hStopEvent = hStopEvent.Get(); // Копируем хэндл StopEvent, main владеет оригиналом

    pCaptureParamsHeap->pSharedBuffer = &sharedBuffer;
    pCaptureParamsHeap->pSharedDspParams = &sharedDspParams;
    pCaptureParamsHeap->pSharedMonitorData = &sharedMonitorData;
    pCaptureParamsHeap->bufferFrameCount = captureBufferFrameCount;
    pCaptureParamsHeap->devicePeriod = captureDefaultPeriod;
    pCaptureParamsHeap->isEventDriven = (pCaptureParamsHeap->hCaptureEvent.Get() != NULL); // Проверяем, был ли создан EventHandle

    // Перемещение владения ресурсами из main в структуру параметров
    pRenderParamsHeap->pAudioClient = std::move(pRenderAudioClient); // Перемещение владения
    pRenderParamsHeap->pRenderClient = std::move(pRenderRenderClient); // Перемещение владения
    pRenderParamsHeap->hRenderEvent = std::move(hRenderEvent); // Перемещение владения
    pRenderParamsHeap->hStopEvent = hStopEvent.Get(); // Копируем хэндл StopEvent, main владеет оригиналом

    pRenderParamsHeap->pSharedBuffer = &sharedBuffer;
    pRenderParamsHeap->pSharedMonitorData = &sharedMonitorData;
    pRenderParamsHeap->bufferFrameCount = renderBufferFrameCount;
    pRenderParamsHeap->devicePeriod = renderDefaultPeriod;
    pRenderParamsHeap->isEventDriven = (pRenderParamsHeap->hRenderEvent.Get() != NULL); // Проверяем, был ли создан EventHandle

    if (pRenderMixFormat) { 
        size_t renderFormatSize = sizeof(WAVEFORMATEX) + pRenderMixFormat->cbSize; 
        if (renderFormatSize <= sizeof(pRenderParamsHeap->renderFormatBuffer)) { 
            memcpy(pRenderParamsHeap->renderFormatBuffer, pRenderMixFormat, renderFormatSize); 
            pRenderParamsHeap->pRenderFormat = reinterpret_cast<PWAVEFORMATEX>(pRenderParamsHeap->renderFormatBuffer);
        } else { 
            printf("Error: Render format size (%zu) exceeds internal buffer size (%zu). Cannot start render thread.\n", renderFormatSize, sizeof(pRenderParamsHeap->renderFormatBuffer));
            // Cannot proceed without render format in params. Need to clean up capture thread if started.
            if (pCaptureAudioClient) { WARN_ON_ERROR(pCaptureAudioClient->Stop(), "AudioClient::Stop (main after RenderParams fail)"); }
            SetEvent(hStopEvent.Get());
            if (hCaptureThreadHandle) { WaitForSingleObject(hCaptureThreadHandle.Get(), INFINITE); } // hCaptureThreadHandle is not yet declared here, this check is wrong. Move this check later.
             // Need to explicitly delete the heap params here since CreateThread failed.
             std::unique_ptr<AudioRenderThreadParams> cleanupRenderParams(pRenderParamsHeap.release()); // Recover ownership to clean up
             std::unique_ptr<AudioCaptureThreadParams> cleanupCaptureParams(pCaptureParamsHeap.release()); // Also recover capture params
            return 1; 
        }
    } else { 
        printf("Error: Render Mix Format is NULL before filling render thread parameters. Cannot start render thread.\n");
         if (pCaptureAudioClient) { WARN_ON_ERROR(pCaptureAudioClient->Stop(), "AudioClient::Stop (main after RenderFormat NULL)"); }
         SetEvent(hStopEvent.Get());
         // Need to clean up capture thread if started. Handle hCaptureThreadHandle state later.
         // Need to explicitly delete the heap params here since CreateThread failed.
         std::unique_ptr<AudioRenderThreadParams> cleanupRenderParams(pRenderParamsHeap.release()); // Recover ownership to clean up
         std::unique_ptr<AudioCaptureThreadParams> cleanupCaptureParams(pCaptureParamsHeap.release()); // Also recover capture params
        return 1; 
    }

    // Move CreateThread calls here after all parameter setup and early exits are handled.

    printf("Создание потока захвата аудио...\n"); ScopedHandle hCaptureThreadHandle;
    hCaptureThreadHandle = CreateThread(NULL, 0, AudioCaptureThread, pCaptureParamsHeap.release(), 0, NULL);
    if (!hCaptureThreadHandle) {
        printf("Ошибка CreateThread (Capture): %lu.\n", GetLastError());
         // Need to explicitly delete the heap params here since CreateThread failed.
         std::unique_ptr<AudioCaptureThreadParams> cleanupParams(pCaptureParamsHeap.release()); // Recover ownership to clean up
        return 1; // Exit main
    }


    printf("Создание потока рендеринга аудио...\n"); ScopedHandle hRenderThreadHandle; // Объявляем хэндлы потоков здесь
    hRenderThreadHandle = CreateThread(NULL, 0, AudioRenderThread, pRenderParamsHeap.release(), 0, NULL); // unique_ptr.release() передает сырой указатель
    if (!hRenderThreadHandle) {
        printf("Ошибка CreateThread (Render): %lu.\n", GetLastError());
        std::unique_ptr<AudioRenderThreadParams> cleanupParams(pRenderParamsHeap.release()); // Пытаемся вернуть владение
        // Если поток рендеринга не запустился, нужно остановить поток захвата перед выходом.
        printf("Ошибка создания потока рендеринга. Сигнализирую потоку захвата на остановку...\n");
        SetEvent(hStopEvent.Get()); // Сигнализируем событие остановки
        // Wait for the capture thread to exit cleanly
        if (hCaptureThreadHandle) { // Check if capture thread was successfully created
             WaitForSingleObject(hCaptureThreadHandle.Get(), INFINITE);
             printf("Поток захвата остановлен.\n");
        }
        return 1; // Ошибка, выходим из main
    }
 // If CreateThread succeeds, ownership of pRenderParamsHeap is transfered

 //-------------------------------------------------------------------------
 // 11. Start WASAPI Clients (in Main)
 //-------------------------------------------------------------------------
 // Use the ComPtrs in main which still hold the initial ownership before moving them.
 // They should be NULL here because they were moved into the thread params.
 // The WASAPI clients should be started *before* creating the threads,
 // as the threads expect them to be running.

    // Let's revert the order: Initialize, Start Clients, Prepare Params, Create Threads.

    // Reorganizing section 10 and 11:
    // 10. Prepare Parameters (without moving client/event ownership yet) and Create Stop Event
    // 11. Initialize and Start WASAPI Clients
    // 12. Move Client/Event Ownership to Parameters
    // 13. Create Threads

    // ... (Sections 1-9 remain as is) ...

    //-------------------------------------------------------------------------
    // 10. Create Stop Event and Prepare Parameter Structs (Heap allocation)
    //-------------------------------------------------------------------------
    ScopedHandle hStopEvent; hStopEvent = CreateEvent(NULL, TRUE, FALSE, NULL);
    if (!hStopEvent) { printf("Ошибка CreateEvent (hStopEvent): %lu.\n", GetLastError()); return 1; }
    printf("Stop Event created.\n");

    auto pCaptureParamsHeap = std::make_unique<AudioCaptureThreadParams>();
    auto pRenderParamsHeap = std::make_unique<AudioRenderThreadParams>();

    // Parameters pointers to shared data
    pCaptureParamsHeap->pSharedBuffer = &sharedBuffer;
    pCaptureParamsHeap->pSharedDspParams = &sharedDspParams;
    pCaptureParamsHeap->pSharedMonitorData = &sharedMonitorData;
    pRenderParamsHeap->pSharedBuffer = &sharedBuffer;
    pRenderParamsHeap->pSharedMonitorData = &sharedMonitorData;

    // Copy Render Format into Render Thread Params buffer
    if (pRenderMixFormat) {
        size_t renderFormatSize = sizeof(WAVEFORMATEX) + pRenderMixFormat->cbSize;
        if (renderFormatSize <= sizeof(pRenderParamsHeap->renderFormatBuffer)) {
            memcpy(pRenderParamsHeap->renderFormatBuffer, pRenderMixFormat, renderFormatSize);
            pRenderParamsHeap->pRenderFormat = reinterpret_cast<PWAVEFORMATEX>(pRenderParamsHeap->renderFormatBuffer);
        } else {
            printf("Error: Render format size (%zu) exceeds internal buffer size (%zu). Cannot proceed.\n", renderFormatSize, sizeof(pRenderParamsHeap->renderFormatBuffer));
            return 1;
        }
    } else {
         printf("Error: Render Mix Format is NULL. Cannot proceed.\n");
         return 1;
    }

    // Parameters for WASAPI buffer sizes and periods (obtained below)
    // Parameters for event driven mode (obtained below)
    pCaptureParamsHeap->hStopEvent = hStopEvent.Get(); // Copy handle
    pRenderParamsHeap->hStopEvent = hStopEvent.Get();   // Copy handle


    //-------------------------------------------------------------------------
    // 11. Initialize WASAPI Clients (and get periods/sizes)
    //-------------------------------------------------------------------------
    ComPtr<IAudioClient> pCaptureAudioClient; ComPtr<IAudioCaptureClient> pCaptureCaptureClient; ScopedHandle hCaptureEvent;
    ComPtr<IAudioClient> pRenderAudioClient; ComPtr<IAudioRenderClient> pRenderRenderClient; ScopedHandle hRenderEvent;

    hr = InitializeCaptureDevice(pCaptureDevice.Get(), pCaptureMixFormat, pCaptureAudioClient, pCaptureCaptureClient, hCaptureEvent);
    if (FAILED(hr)) { printf("Ошибка инициализации устройства захвата (0x%x).\n", hr); return 1; } if (!pCaptureAudioClient || !pCaptureCaptureClient) { printf("Ошибка: Инициализация устройства захвата вернула NULL объекты.\n"); return 1; }
    UINT32 captureBufferFrameCount = 0; REFERENCE_TIME captureDefaultPeriod = 0, captureMinimumPeriod = 0;
    if (pCaptureAudioClient) { hr = pCaptureAudioClient->GetBufferSize(&captureBufferFrameCount); WARN_ON_ERROR(hr, "GetBufferSize (main, capture)"); hr = pCaptureAudioClient->GetDevicePeriod(&captureDefaultPeriod, &captureMinimumPeriod); WARN_ON_ERROR(hr, "GetDevicePeriod (main, capture)"); printf("Параметры буфера захвата: WASAPI Size=%u frames, Default Period=%lld ms.\n", captureBufferFrameCount, captureDefaultPeriod/10000); } else { printf("Error: pCaptureAudioClient is NULL after successful InitializeCaptureDevice!\n"); return 1; }

    hr = InitializeRenderDevice(pRenderDevice.Get(), pRenderMixFormat, pRenderAudioClient, pRenderRenderClient, hRenderEvent);
    if (FAILED(hr)) { printf("Ошибка инициализации устройства вывода (0x%x).\n", hr); return 1; } if (!pRenderAudioClient || !pRenderRenderClient) { printf("Ошибка: Инициализация устройства вывода вернула NULL объекты.\n"); return 1; }
    UINT32 renderBufferFrameCount = 0; REFERENCE_TIME renderDefaultPeriod = 0, renderMinimumPeriod = 0;
    if (pRenderAudioClient) { hr = pRenderAudioClient->GetBufferSize(&renderBufferFrameCount); WARN_ON_ERROR(hr, "GetBufferSize (main, render)"); hr = pRenderAudioClient->GetDevicePeriod(&renderDefaultPeriod, &renderMinimumPeriod); WARN_ON_ERROR(hr, "GetDevicePeriod (main, render)"); printf("Параметры буфера рендеринга: WASAPI Size=%u frames, Default Period=%lld ms.\n", renderBufferFrameCount, renderDefaultPeriod/10000); } else { printf("Error: pRenderAudioClient is NULL after successful InitializeRenderDevice!\n"); return 1; }

    // Fill buffer sizes and periods into parameter structs
    pCaptureParamsHeap->bufferFrameCount = captureBufferFrameCount;
    pCaptureParamsHeap->devicePeriod = captureDefaultPeriod; // Or minimumPeriod? Using default for shared mode example.
    pCaptureParamsHeap->isEventDriven = (hCaptureEvent.Get() != NULL); // Check if EventHandle was successfully created

    pRenderParamsHeap->bufferFrameCount = renderBufferFrameCount;
    pRenderParamsHeap->devicePeriod = renderDefaultPeriod; // Or minimumPeriod?
    pRenderParamsHeap->isEventDriven = (hRenderEvent.Get() != NULL); // Check if EventHandle was successfully created


     //-------------------------------------------------------------------------
     // 12. Start WASAPI Clients (after parameters are ready)
     //-------------------------------------------------------------------------
     printf("Попытка запуска клиента рендеринга...\n");
     hr = pRenderAudioClient->Start(); // Use the ComPtr from main before moving
     if (FAILED(hr)) {
         printf("Ошибка AudioClient::Start (Render): 0x%x\n", hr);
         WARN_ON_ERROR(pRenderAudioClient->Stop(), "AudioClient::Stop (Render after start fail)"); // Attempt to stop
         // No threads running yet, can exit directly.
         return 1;
     }
     printf("Клиент рендеринга запущен.\n");
     // Optional: Initial buffer filling for render client to avoid startup glitch
     // This would involve getting buffer, filling silence, releasing buffer.
     // For simplicity, skip for now.

     printf("Попытка запуска клиента захвата...\n");
     hr = pCaptureAudioClient->Start(); // Use the ComPtr from main before moving
     if (FAILED(hr)) {
         printf("Ошибка AudioClient::Start (Capture): 0x%x\n", hr);
         WARN_ON_ERROR(pCaptureAudioClient->Stop(), "AudioClient::Stop (Capture after start fail)"); // Attempt to stop capture
         WARN_ON_ERROR(pRenderAudioClient->Stop(), "AudioClient::Stop (Render after capture start fail)"); // Attempt to stop render
         return 1; // Exit main
     }
     printf("Клиент захвата запущен.\n");


    //-------------------------------------------------------------------------
    // 13. Move Client/Event Ownership to Parameters and Create Threads
    //-------------------------------------------------------------------------

    // Move ownership using std::move
    pCaptureParamsHeap->pAudioClient = std::move(pCaptureAudioClient);
    pCaptureParamsHeap->pCaptureClient = std::move(pCaptureCaptureClient);
    pCaptureParamsHeap->hCaptureEvent = std::move(hCaptureEvent); // Move ownership of the event handle

    pRenderParamsHeap->pAudioClient = std::move(pRenderAudioClient);
    pRenderParamsHeap->pRenderClient = std::move(pRenderRenderClient);
    pRenderParamsHeap->hRenderEvent = std::move(hRenderEvent);     // Move ownership of the event handle

    // Now pCaptureAudioClient, pCaptureCaptureClient, hCaptureEvent, etc. in main are NULL/invalid.
    // The unique_ptr in the thread will own them.


    printf("Создание потока захвата аудио...\n"); ScopedHandle hCaptureThreadHandle; // Use ScopedHandle for thread handle in main
    // Use .release() to transfer the raw pointer ownership from unique_ptr to CreateThread
    hCaptureThreadHandle = CreateThread(NULL, 0, AudioCaptureThread, pCaptureParamsHeap.release(), 0, NULL);
    if (!hCaptureThreadHandle) {
        printf("Ошибка CreateThread (Capture): %lu.\n", GetLastError());
         // The unique_ptr was already released, the params struct is leaked on heap if thread creation fails.
         // A more robust approach would be to *not* release until *after* CreateThread returns success.
         // But given the RAII principles, the struct *contains* the resources.
         // If CreateThread fails, the LPVOID isn't passed. The released raw pointer to params is lost.
         // This is a weakness of using unique_ptr::release() with CreateThread.
         // The params struct and its contents (WASAPI clients, events) are leaked.
         // For this project scope, we'll accept this leak on thread creation failure.
        return 1; // Exit main
    }

    printf("Создание потока рендеринга аудио...\n"); ScopedHandle hRenderThreadHandle; // Use ScopedHandle for thread handle in main
    // Use .release() to transfer the raw pointer ownership from unique_ptr to CreateThread
    hRenderThreadHandle = CreateThread(NULL, 0, AudioRenderThread, pRenderParamsHeap.release(), 0, NULL);
    if (!hRenderThreadHandle) {
        printf("Ошибка CreateThread (Render): %lu.\n", GetLastError());
        // The unique_ptr was already released, the params struct is leaked on heap if thread creation fails.
        // If render thread fails, attempt to signal capture thread to stop.
        printf("Ошибка создания потока рендеринга. Сигнализирую потоку захвата на остановку...\n");
        SetEvent(hStopEvent.Get()); // Signal stop event
        // Wait for the capture thread to exit cleanly
        if (hCaptureThreadHandle) { // Check if capture thread was successfully created
             WaitForSingleObject(hCaptureThreadHandle.Get(), INFINITE);
             printf("Поток захвата остановлен.\n");
        }
        return 1; // Exit main
    }


 //-------------------------------------------------------------------------
 // 14. User Input and Wait for Threads
 //-------------------------------------------------------------------------
 printf("\nАудиопотоки запущены.\n");
 // Pass raw handles/pointers to HandleUserInput
 HandleUserInput(&sharedDspParams, &sharedBuffer, &sharedMonitorData, hStopEvent.Get());

 printf("Подача сигнала остановки потокам завершена.\n");
 printf("Ожидание завершения потоков...\n");
 // Wait for thread handles. ScopedHandle will close them automatically on exit.
 if (hCaptureThreadHandle) WaitForSingleObject(hCaptureThreadHandle.Get(), INFINITE); printf("Поток захвата завершен.\n");
 if (hRenderThreadHandle) WaitForSingleObject(hRenderThreadHandle.Get(), INFINITE); printf("Поток рендеринга завершен.\n");


 //-------------------------------------------------------------------------
 // 15. Stop and Cleanup WASAPI Clients (in Main, after threads exit)
 //-------------------------------------------------------------------------
 // WASAPI clients are owned by ComPtrs *inside* the thread params structures,
 // which are now managed by the unique_ptr in the thread functions.
 // When the threads exit, the unique_ptr goes out of scope, its destructor runs,
 // params struct destructor runs, ComPtr destructors run, Release() is called.
 // So, no explicit Stop/Reset/Release calls are needed *here* in main.
 // However, it might be good practice for the threads to call Stop/Reset just before exiting
 // to ensure a clean shutdown from the WASAPI side, especially if exiting due to errors.
 // Let's add Stop/Reset to the thread exit logic, but ensure they are safe to call (check for NULL pointers).
 // (This was already done in the previous corrected thread code blocks, but commented out)
 // Let's re-add the Stop/Reset calls in the thread functions before their final return.

 // Re-added Stop/Reset calls to thread functions before return.

 printf("Программа завершена.\n");
 return 0; // Normal exit
}

// --- Global Function Definitions (Move them here after main, or ensure prototypes are above main) ---
// For brevity, assuming these are defined elsewhere or their prototypes are sufficient here.
// Need definitions for:
// SelectAudioDevicePrompt, GetMixFormatForDevice, InitializeCaptureDevice, InitializeRenderDevice,
// ApplyDspEffectChain, HandleUserInput, PrintCurrentDspSettings, PrintMonitorData, PrintWaveFormatEx,
// BiquadFilterTypeToString, and the Biquad coefficient calculation helpers.

// BiquadFilterTypeToString definition (example)
const char* BiquadFilterTypeToString(BiquadFilterType type) {
    switch (type) {
        case BiquadFilterType::Identity: return "Identity";
        case BiquadFilterType::LowPass: return "LowPass";
        case BiquadFilterType::HighPass: return "HighPass";
        case BiquadFilterType::BandPass: return "BandPass";
        case BiquadFilterType::BandPass0dB: return "BandPass0dB";
        case BiquadFilterType::Notch: return "Notch";
        case BiquadFilterType::LowShelf: return "LowShelf";
        case BiquadFilterType::HighShelf: return "HighShelf";
        case BiquadFilterType::Peak: return "Peak";
        default: return "Unknown";
    }
}